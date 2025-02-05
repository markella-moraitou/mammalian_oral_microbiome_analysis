##### IDENTIFY CORE MICROBIOMES #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(ggVennDiagram)
library(microViz)
library(vegan)
library(rphylopic)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "core_microbiomes")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

source(file.path("ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#################
#### PROCESS ####
#################

## Calculate prevalence of each OTU in each species

# Prevalence per species
species <- phy_sp_f@sam_data$Species %>% levels
prevalence <- data.frame(row.names = species)

for (spe in species) {
  subset <- phy_sp_f %>% subset_samples(Species == spe)
  temp <- prevalence(subset, detection = 2/100, , include.lowest = TRUE, sort = TRUE) %>% data.frame() %>% t
  rownames(temp) <- spe
  prevalence <- rbind(prevalence, temp)
}

prevalence <- t(prevalence) %>% data.frame %>% arrange(desc(rowSums(.))) %>%
                rownames_to_column("taxon")

write.csv(prevalence, file = file.path(subdir, "taxa_prevalence.csv"), quote = FALSE, row.names = FALSE)

#### Get core taxa per species ####
## 90% prevalence

core_species <- data.frame(core_taxa = character(), host_species = character())

for (sp in species) {
  sp_dot <- gsub(pattern = " ", replacement = ".", sp)
  # Get core microbiota for this order
  core_taxa <- prevalence %>% select(taxon, !!sym(sp_dot)) %>%
        filter(!!sym(sp_dot) >= 0.9) %>% pull(taxon)
  df <- data.frame(core_taxa = core_taxa,
                   host_species = sp)
  core_species <- rbind(core_species, df)
}

# Add extra information
# on the taxa
core_species$core_genus <- as.vector(phy_sp_f@tax_table[match(core_species$core_taxa, taxa_names(phy_sp_f)), "genus"])

# and on the samples
core_species$host_order <- phy_sp_f@sam_data$Order[match(core_species$host_species, phy_sp_f@sam_data$Species)]

core_species <- core_species %>% arrange(host_order, host_species, core_taxa) %>% 
        select(host_order, host_species, core_genus, core_taxa)

# Calculate number of core taxa per species
core_species_summary <- 
        core_species %>% group_by(host_species, host_order) %>% 
        summarise(n_core_taxa = n(),
               n_core_genus = n_distinct(core_genus))

# Save tables
write.csv(core_species, file = file.path(subdir, "core_mb_per_host_species.csv"), quote = FALSE, row.names = FALSE)
write.csv(core_species_summary, file = file.path(subdir, "core_mb_per_host_species_summary.csv"), quote = FALSE, row.names = FALSE)

# For orders with more than two species, get core taxa per order
# defined as taxa that are considered core in at least 80% of species in that order

core_species_per_order <- core_species %>% group_by(host_order) %>%
        filter(n_distinct(host_species) > 2) %>%
        # Calculate number of species in order
        mutate(n_species = n_distinct(host_species)) %>%
        # Calculate prevalence in order
        group_by(core_taxa, core_genus, host_order) %>%
        summarise(prevalence_in_order = n_distinct(host_species)/first(n_species)) %>%
        filter(prevalence_in_order >= 0.9) %>%
        arrange(host_order, core_taxa) %>% ungroup

# Calculate number of core taxa per order
core_species_per_order_summary <- core_species_per_order %>% group_by(host_order) %>%
        summarise(n_core_taxa = n(),
               n_core_genus = n_distinct(core_genus))

# Save tables
write.csv(core_species_per_order, file = file.path(subdir, "core_mb_per_host_order.csv"), quote = FALSE, row.names = FALSE)
write.csv(core_species_per_order_summary, file = file.path(subdir, "core_mb_per_host_order_summary.csv"), quote = FALSE, row.names = FALSE)

#### Plot Venn Diagram of genera ####

# Collect core genera in a list
core_genera <- list()

for (ord in unique(core_species_per_order$host_order)) {
    core_genera[[ord]] <-
        core_species_per_order %>% filter(host_order == ord) %>%
        pull(core_genus) %>% unique
}

core_genera_venn <- 
  ggVennDiagram(core_genera, set_color = order_palette[names(core_genera)], label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. genera") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom")

ggsave(core_genera_venn, file=file.path(subdir, "core_genera_venn.png"),
       device = "png", width = 8, height = 8)

#### PCA with only core order taxa ####

phy_core <- phy_sp_f %>% subset_taxa(taxa_names(phy_sp_f) %in% core_species_per_order$core_taxa) %>%
    subset_samples(Order %in% core_species_per_order$host_order) %>%
    # Agglomerate by genus
    tax_glom(taxrank = "genus")

taxa_names(phy_core) <- phy_core@tax_table[,"genus"]

phy_core_clr <- transform(phy_core, "clr")

# PCA ordination
ord <- ord_calc(phy_core_clr, method = "PCA") 

# Get shape scales for plotting
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

# Axes 1 & 2
p <- 
  ord_plot(ord, colour="Order", shape = "diet.general", alpha = 0.5, plot_taxa = 1:8) +
  custom_theme() +
  scale_shape_manual(values = diet_shape_scale, name = "Diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_core_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_core_clr)$uid, width = 0.07, alpha = 0.8)

ggsave(p, file = file.path(subdir, "pca_core_order_1_2.png"), width = 8, height = 8)

# Axes 3 & 4
p <- 
  ord_plot(ord, colour="Order", shape = "diet.general", alpha = 0.5, plot_taxa = 1:8, axes = c(3, 4)) +
  custom_theme() +
  scale_color_manual(values=order_palette, name = "Order") +
  scale_shape_manual(values = diet_shape_scale, name = "Diet") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_core_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_core_clr)$uid, width = 0.07, alpha = 0.8)

ggsave(p, file = file.path(subdir, "pca_core_order_3_4.png"), width = 8, height = 8)

#### Heatmap ####

# Get more info on taxon category
heat_data <- transform(phy_core, "log10p") %>% psmelt %>%
        select(OTU, Sample, Order, Abundance) %>%
        mutate(taxon_type = case_when(OTU %in% Reduce(intersect, core_genera) ~ "Mammalian core",
                                      OTU %in% setdiff(core_genera[["Carnivora"]], unlist(core_genera[!names(core_genera) %in% "Carnivora"])) ~ "Carnivora core only",
                                      OTU %in% intersect(core_genera[["Carnivora"]], intersect(core_genera[["Perissodactyla"]], core_genera[["Primates"]])) ~ "Non artiodactyl core",
                                      TRUE ~ "Other")) %>%
        mutate(taxon_type = factor(taxon_type, levels = c("Mammalian core", "Non artiodactyl core", "Carnivora core only", "Other"))) %>%
        # Shorten order names for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

p <- ggplot(heat_data, aes(x = Sample, y = OTU, fill = Abundance)) +
        geom_tile() +
        facet_grid(cols = vars(Order), scales = "free", space = "free",
                   rows = vars(taxon_type), switch = "y") +
        scale_fill_gradientn(colors = grad_palette, name = "Abundance (log 10)") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "bottom",
              panel.spacing.y = unit(1.5, "lines"),
              strip.placement = "outside",
              strip.clip = "off",
              strip.text.y.left = element_text(angle=0, vjust=1),
              strip.text.y = element_text(margin = margin(t=-15, r = -110),
                                          size = 15, hjust = 1),
              strip.background.y = element_blank())

ggsave(file.path(subdir, "heatmap_core_genera.png"), p, width=8, height=9)

#### Percentage of abundance consisting of core taxa ####
# This reflects how representative the above PCA is of the whole community

# Calculate the proportion of abundance that is made up of core taxa
phy_melt <- transform(phy_sp_f, "compositional") %>% psmelt %>% filter(Abundance > 0)

core_species_per_order$is.core <- "Core in this order"

core_taxa_abund <-
        phy_melt %>%
        select(OTU, genus, Sample, Species, Common.name, Order, Abundance) %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_species_per_order$host_order) %>%
        # Get core taxa per order
        left_join(core_species_per_order, by = c("OTU" = "core_taxa", "Order" = "host_order", "genus" = "core_genus"),
                  relationship = "many-to-one") %>%
        # If a taxon isnt core in that order, check if it is core in another order
        mutate(is.core = case_when(!is.na(is.core) ~ is.core,
                                   OTU %in% core_species_per_order$core_taxa ~ "Core in other order")) %>%
        mutate(is.core = factor(is.core, levels = c("Core in other order", "Core in this order"))) %>%
        filter(!is.na(is.core)) %>%
        # Get relative abunddance of host taxa per sample
        group_by(Sample, Species, Common.name, Order, is.core) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Summarize by species
        group_by(Species, Common.name, Order, is.core) %>%
        summarise(mean_abundance = mean(Abundance),
                  sd = sd(Abundance, na.rm = TRUE)) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Plot
p <- ggplot(core_taxa_abund, aes(y = Common.name, x = mean_abundance, fill = is.core)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#DA8A3D", "#DA4C3D"), name = "") +
        facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
        labs(x = "Relative abundance of core taxa") +
        theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = guide_legend(reverse = TRUE))

ggsave(p, file = file.path(subdir, "core_taxa_abundance.png"), width = 8, height = 8)

