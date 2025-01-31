##### EXPLORE FILTERED & DECONTAMINATED DATA #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(microbiomeutilities)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "filtered_data")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

#########################
#####  COMPOSITION  #####
#########################

#### PHYLUM COMPOSITION ####
phy_phylum <- tax_glom(phy_sp_f, taxrank = "phylum")
taxa_names(phy_phylum) <- as.vector(phy_phylum@tax_table[,"phylum"])

# Get only 10 most common phyla and turn rest to other
phylum_grouped = data.frame(abundance = taxa_sums(phy_phylum), superkingdom = phy_phylum@tax_table[,"superkingdom"]) %>% rownames_to_column("phylum") %>% arrange(-abundance) %>%
  mutate(phylum_grouped = ifelse(row_number() > 5, paste("Other", superkingdom, sep = " "), phylum))

# Change phylum names to grouped names and reaggregate
phy_phylum@tax_table[,"phylum"] <- phylum_grouped$phylum_grouped[match(phy_phylum@tax_table[,"phylum"], phylum_grouped$phylum)]

# Aggregate again
phy_phylum <- tax_glom(phy_phylum, taxrank = "phylum")
taxa_names(phy_phylum) <- phy_phylum@tax_table[,"phylum"] 

# Melt and turn phyla into a factor and reorder
phy_phylum_melt <- psmelt(transform(phy_phylum, "compositional"))
phy_phylum_melt$OTU <- factor(phy_phylum_melt$OTU , levels=names(phylum_palette))

# Order by fusobacteriota
sample_levels <- select(phy_phylum_melt, c(Sample, Species, Order_grouped, OTU, Abundance)) %>% filter(OTU == "Pseudomonadota") %>%
  arrange(Order_grouped, Species, desc(Abundance)) %>% pull(Sample)

phy_phylum_melt$Sample <- factor(phy_phylum_melt$Sample, levels=sample_levels)

p = ggplot(data = phy_phylum_melt, aes(x = Abundance, y = Sample, fill = OTU)) +
  geom_bar(stat = "identity") +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  scale_fill_manual(values=phylum_palette, name = "Phylum") +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title.position = "top", legend.key.spacing.x = unit(0.5, "cm"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  xlab("")

## Get sample_metadata
species_bar <- phy_phylum_melt %>% select(Sample, Species, Common.name, Order_grouped) %>% unique %>%
               arrange(Sample) %>%
               # Choose one label per species (in the middle of the species group)
               group_by(Species) %>% mutate(order = row_number()) %>%
               mutate(label = ifelse(order == floor(mean(order)), Common.name, "")) %>%
               mutate(label = ifelse(is.na(label), as.character(Species), label))

p_bar <-
  ggplot(data = species_bar, aes(y = Sample, x=1, fill = Species, group = Common.name)) +
  geom_tile() + scale_fill_manual(values = species_palette, name = "") +
  facet_grid(rows = vars(Order_grouped), scales = "free", space = "free", switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank()) + scale_y_discrete(position = "right", label = setNames(species_bar$label, species_bar$Sample)) +
    xlab("") + ylab("")

ggsave(filename = file.path(subdir, "phy_sp_f_composition.png"), device="png", width=12, height=20,
       plot_grid(p, p_bar, ncol = 2, align = "h", rel_widths = c(4.5, 1.5)))

####################
#### ORDINATION ####
####################

#### All data ####
ord <- ordinate(phy_sp_f_norm, "PCoA", "euclidean")

# Get shape scales for plotting
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)
species_shape_scale <- 1:length(unique(phy_sp_f_norm@sam_data$Common.name))
names(species_shape_scale) <- unique(phy_sp_f_norm@sam_data$Common.name)

# Color by order
p <- plot_ordination(phy_sp_f_norm, ord, color="Order_grouped", shape="calculated_species_main_diet") +
  geom_point(size=2, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_order_1_2.png"), p, width=8, height=6)

p <- plot_ordination(phy_sp_f_norm, ord, axes = c(3, 4), color="Order_grouped", shape="calculated_species_main_diet") +
  geom_point(size=2, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_order_3_4.png"), p, width=8, height=6)

# Color by diet
p <- plot_ordination(phy_sp_f_norm, ord, color="calculated_species_main_diet", shape="Order_grouped") +
  geom_point(size=2, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_diet_1_2.png"), p, width=8, height=6)

p <- plot_ordination(phy_sp_f_norm, ord, axes = c(3, 4), color="calculated_species_main_diet", shape="Order_grouped") +
  geom_point(size=2, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_diet_3_4.png"), p, width=8, height=6)

#### Deep dataset ####
ord <- ordinate(phy_deep_norm, "PCoA", "euclidean")

p <- plot_ordination(phy_deep_norm, ord, color="Order", shape="calculated_species_main_diet") +
  geom_point(size=2, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_subset_deep_1_2.png"), p, width=8, height=6)

#### Artiodactyla ####
ord <- ordinate(phy_artio_norm, "PCoA", "euclidean")

p <- plot_ordination(phy_artio_norm, ord, color="calculated_species_main_diet", shape="Common.name") +
  geom_point(size=2, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = species_shape_scale, name = "Species") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_subset_artio_1_2.png"), p, width=8, height=6)

#### Carnivora ####
ord <- ordinate(phy_carni_norm, "PCoA", "euclidean")

p <- plot_ordination(phy_carni_norm, ord, color="calculated_species_main_diet", shape="Common.name") +
  geom_point(size=3, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = c(15, 2, 4)) +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_subset_carni_1_2.png"), p, width=8, height=6)

#### Primates ####
ord <- ordinate(phy_prim_norm, "PCoA", "euclidean")

p <- plot_ordination(phy_prim_norm, ord, color="calculated_species_main_diet", shape="Common.name") +
  geom_point(size=3, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = 1:length(levels(phy_prim_norm@sam_data$Species))) +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_subset_prim_1_2.png"), p, width=8, height=6)

#### Habitat ####
ord <- ordinate(phy_habitat_norm, "PCoA", "euclidean")

p <- plot_ordination(phy_habitat_norm, ord, color="habitat.general", shape = "Common.name") +
  geom_point(size=3, alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = species_shape_scale) +
  scale_color_manual(values = habitat_palette, name = "Habitat") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "PCoA_subset_habitat_1_2.png"), p, width=8, height=6)

#################
#### HEATMAP ####
#################

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

png(filename = file.path(subdir, "heatmap_all_taxa.png"), width=16, height=20, units="in", res=300)
plot_taxa_heatmap(phy_sp_f, subset.top=ntaxa(phy_sp_f), transformation="clr",
                  VariableA=c("Order", "calculated_species_main_diet", "unmapped_count"),
                  annotation_colors = list("Order" = order_palette, "calculated_species_main_diet" = diet_palette),
                  show_rownames = FALSE,
                  show_colnames = FALSE, heatcolors = grad_palette)$plot
dev.off()

#########################
#### ALPHA DIVERSITY ####
#########################

# Rarefy to min depth
phy_sp_rarefied <- rarefy_even_depth(subset_samples(phy_sp_f, sample_sums(phy_sp_f) > 100), sample.size = 100, rngseed = 1)

# Create a named vector for the relabeling
order_labels <- c("Rodentia" = "Rod.", "Carnivora" = "Carniv.")

# Plot
p <- plot_richness(phy_sp_f, x="Species", measures=c("Observed")) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp_f@sam_data$Common.name, phy_sp_f@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness") +
  coord_flip()
ggsave(file.path(subdir, "alpha_diversity.png"), p, width=8, height=10)

# Plot
p <- plot_richness(phy_sp_rarefied, x="Species", measures=c("Observed")) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp_f@sam_data$Common.name, phy_sp_f@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness (after rarefaction)") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_rarefied.png"), p, width=8, height=10)

#### RAREFACTION CURVES ####
# get it separate per tax. order
subsamples <- seq(0, 1000, by=100)[-1]
plot_list <- list()
for (ord in unique(phy_sp_f@sam_data$Order_grouped)) { 
  phy_subset <- phy_sp_f %>% subset_samples(Order_grouped == ord)
  p <- plot_alpha_rcurve(phy_subset, index="observed", subsamples=subsamples,
                  lower.conf = 0.025, upper.conf = 0.975,
                  group="Species", label.color = "brown3",
                  label.size = 3, label.min = TRUE) +
        scale_fill_manual(values=species_palette, name = "Species") +
        scale_color_manual(values=species_palette, name = "Species")
  plot_list[[ord]] <- p
}

p <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

ggsave(file.path(subdir, "rarefaction_curves.png"), p, width=8, height=16)
