##### PLOT SUMMARY PLOTS FOR PHYLOSEQ OBJECTS #####

#### Various plots including abundance-prevalence, composition, rank abundance, taxonomic richness, ordinations and heatmap

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)
library(microViz)
library(microbiomeutilities)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "raw_data_plots")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

source("ordination_functions.R")

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp <- readRDS(file.path(phydir, "phy_sp.RDS"))
phy_sp_clr <- readRDS(file.path(phydir, "phy_sp_clr.RDS"))

#######################
#### SUMMARY PLOTS ####
#######################

phy_spr <- transform(phy_sp, "compositional")

#### Composition plot ####
phy_phylum <- tax_glom(phy_sp, taxrank = "phylum")
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

## Get info on sequencing depth
rc <- phy_phylum_melt %>% select(Sample, Species, Order_grouped, unmapped_count) %>% arrange(Sample)

p_count <- ggplot(data = rc, aes(x = unmapped_count, y = Sample)) +
  geom_point() +
  scale_x_continuous(trans = "log10", breaks = c(0, 10^2, 10^4, 10^6, 10^8)) +
  facet_grid(rows = vars(Order_grouped), scales = "free", space = "free") +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  xlab("Read count")

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

ggsave(filename = file.path(subdir, "phy_sp_composition.png"), device="png", width=12, height=20,
       plot_grid(p, p_count, p_bar, ncol = 3, align = "h", rel_widths = c(4.5, 1, 1.5)))

#### Ordinations ####

# Get sample type as a factor
phy_sp_clr@sam_data$sample_type <- factor(ifelse(!phy_sp_clr@sam_data$is.neg, "sample",
                                             ifelse(phy_sp_clr@sam_data$Species == "Environmental control", "swab", "blank")),
                                             levels = c("sample", "swab", "blank"))

ord <-  ord_calc(phy_sp_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot.png"), p, width=8, height=6)

# Axes 1 and 2
p <- ord_plot(ord, colour="Order_grouped", shape="sample_type", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=c("swab"=0, "blank"=2, "sample"=16), name = "Is control/blank") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left")

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "phy_sp_sample_PCA_1_2.png"), p, width=8, height=10)

# Axes 3 and 4
p <- ord_plot(ord, colour="Order_grouped", shape="sample_type", alpha = 0.5, axes = c(3,4)) +
  custom_theme() +
  scale_shape_manual(values=c("swab"=0, "blank"=2, "sample"=16), name = "Is control/blank") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left")

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(3, 4), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "phy_sp_sample_PCA_3_4.png"), p, width=8, height=10)

#### Heatmap ####
grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

png(filename = file.path(subdir, "phy_sp_heatmap.png"), width=16, height=20, units="in", res=300)
plot_taxa_heatmap(phy_sp, subset.top=ntaxa(phy_sp), transformation="clr",
                  VariableA=c("Order", "diet.general", "unmapped_count"),
                  annotation_colors = list("Order" = order_palette, "diet.general" = diet_palette),
                  show_rownames = FALSE,
                  show_colnames = FALSE, heatcolors = grad_palette)$plot
dev.off()

#### Rarefaction curves ####
#p <- plot_alpha_rcurve(phy_sp) +
#  scale_fill_manual(values = order_palette, name = "Host order") +
#  scale_colour_manual(values = order_palette, name = "Host order")

#ggsave(file.path(subdir, "phy_sp_rarefaction.png"), p, width=8, height=6)

#### Prevalence at different abundance thresholds ####
prevalence <- data.frame()
for (abundance in c(0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) {
  df <- data.frame(prevalence(subset_samples(phy_spr, !is.neg), detection = abundance)) %>% rownames_to_column()
  colnames(df) <- c("OTU", "Prevalence")
  df$Abundance_thres <- abundance
  prevalence <- rbind(prevalence, df)
}

prevalence$Abundance_thres <- as.factor(prevalence$Abundance_thres)

# Get average abundance per OTU
av_abund <- data.frame(OTU = taxa_names(phy_spr), Average_abundance = taxa_sums(phy_spr)/ntaxa(phy_spr)) %>%
  arrange(desc(Average_abundance))

prevalence$OTU <- factor(prevalence$OTU, levels = rev(av_abund$OTU))

# Save table
write.table(prevalence, file.path(subdir, "phy_sp_prevalence.csv"), sep=",", row.names=FALSE, quote=FALSE)

# Plot
p <- ggplot(data = prevalence, aes(x = Abundance_thres, y = OTU, fill = Prevalence)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("Abundance threshold") + ylab("taxa")

ggsave(file.path(subdir, "phy_sp_prevalence.png"), p, width=4, height=8)

#### Sample and taxa sums ####
sample_sums <- data.frame(sum = sample_sums(phy_sp), order = sample_data(phy_sp)$Order_grouped) %>% rownames_to_column("Sample")

ps <- ggplot(data = sample_sums, aes(x = sum, fill = order, colour = order)) +
  geom_density(alpha=0.2, linewidth=1.5) +
  scale_fill_manual(values=order_palette, name = "Order") +
  scale_colour_manual(values=order_palette, name = "Order") +
  scale_x_continuous(trans = "log10") +
  geom_density(data = sample_sums, aes(x = sum), fill = "transparent", linewidth = 1.5, linetype = "dashed", inherit.aes = FALSE)

taxa_sums <- data.frame(sum = taxa_sums(phy_sp)/sum(taxa_sums(phy_sp)), phylum = tax_table(phy_sp)[,"phylum"]) %>% rownames_to_column("Taxon") %>% 
  left_join(select(phylum_grouped, c(phylum, phylum_grouped)), by = "phylum") %>%
  mutate(phylum_grouped = factor(phylum_grouped, levels = names(phylum_palette)))

pt <- ggplot(data = taxa_sums, aes(x = sum, colour = phylum_grouped, fill = phylum_grouped)) +
  geom_density(alpha=0.2, linewidth=1.5) +
  scale_colour_manual(values=phylum_palette, name = "Phylum") +
  scale_fill_manual(values=phylum_palette, name = "Phylum") +
  scale_x_continuous(trans = "log10") +
  geom_density(data = taxa_sums, aes(x = sum), fill = "transparent", linewidth = 1.5, linetype = "dashed", inherit.aes = FALSE)

ggsave(file.path(subdir, "phy_sp_sums_distr.png"),
       plot_grid(ps, pt, ncol = 1, align = "v", rel_heights = c(1, 1)),
       width=8, height=8)

#### Taxonomic richness and read count ####
p = ggplot(phy_sp@sam_data, aes(x=sample_sums(phy_sp), y=taxa_raw, color = is.neg)) +
  geom_point(size=2) +
  labs(x="Number of reads classified", y="Number of OTUs") +
  scale_color_manual(values=c(`TRUE`="grey", `FALSE`="darkgreen"), labels = c(`TRUE`="negative", `FALSE`="sample"), name = "") +
  scale_x_continuous(trans="log10") +
  theme(legend.position = "bottom")

p <- ggMarginal(p, type="histogram", size=2, groupFill=TRUE)

ggsave(file.path(subdir, "phy_sp_taxa_and_read_dist.png"), p, width=4, height=4)

#### Rank abundance plot ####
psm <- phy_spr %>% psmelt

# Get taxon rank within a sample vs abundance for the first 500 taxa
rank_abund <- psm %>% group_by(Sample) %>% mutate(Rank=rank(-Abundance, ties="first")) %>%
  select(OTU, Sample, is.neg, Species, Order, Abundance, Rank, taxa_raw) %>%
  # Keep only the first 500 taxa. Remove 0s
  filter(Rank<500) %>% filter(Abundance>0)

# Save table
write.table(rank_abund, file.path(subdir, "phy_sp_rank_abundance_taxa.csv"), sep=",", row.names=FALSE, quote=FALSE)

# Plot
p <- ggplot(data=rank_abund, aes(x=Rank, y=Abundance, colour=is.neg)) +
  geom_jitter(alpha = 0.3, size = 1, width = 0.05, height = 0.001) +
  scale_y_continuous(trans="log10") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=c(`TRUE`="grey30", `FALSE`="darkgreen"), name = "Is control/blank") +
  ylab("log-transformed relative abundance") + xlab("Abundance rank in sample")

ggsave(file=file.path(subdir, "phy_sp_rank_abundance_plot.png"), p, width=8, height=6)

#### Abundance and prevalence plot ####
set.seed(123)
p <- microbiomeutilities::plot_abund_prev(phy_sp)
ggsave(file.path(subdir, "phy_sp_abund_prev.png"), p, width=8, height=6)
