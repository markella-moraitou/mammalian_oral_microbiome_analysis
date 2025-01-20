##### PLOT SUMMARY PLOTS FOR PHYLOSEQ OBJECTS #####

#### Various plots including abundance-prevalence, composition, rank abundance, taxonomic richness, ordinations and heatmap

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggExtra)
library(microbiomeutilities)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # subdirectory for the output of this script

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp <- readRDS(file.path(subdir, "phyloseq_objects", "phy_sp.RDS"))
phy_sp_norm <- readRDS(file.path(subdir, "phyloseq_objects", "phy_sp_norm.RDS"))

#######################
#### SUMMARY PLOTS ####
#######################

#### Abundance and prevalence plot
set.seed(123)
p <- microbiomeutilities::plot_abund_prev(phy_sp)
ggsave(file.path(subdir, "phy_sp_abund_prev.png"), p, width=8, height=6)

#### Composition plot
phy_phylum <- tax_glom(phy_sp, taxrank = "phylum")
taxa_names(phy_phylum) <- as.vector(phy_phylum@tax_table[,"phylum"])

p = plot_composition(transform(phy_phylum, "compositional"), group_by = "Order") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank())

ggsave(file.path(subdir, "phy_sp_composition.png"), p, width=8, height=6)

#### Rank abundance plot
psm <- phy_sp %>% transform('compositional') %>% psmelt

# Get taxon rank within a sample vs abundance for the first 500 taxa
rank_abund <- psm %>% group_by(Sample) %>% mutate(Rank=rank(-Abundance, ties="first")) %>%
  select(OTU, Sample, is.neg, Species, Order, Abundance, Rank, taxa_raw) %>%
  # Keep only the first 500 taxa. Remove 0s
  filter(Rank<500) %>% filter(Abundance>0)

# Save table
write.table(rank_abund, file.path(subdir, "phy_sp_rank_abundance_taxa.csv"), sep=",", row.names=FALSE, quote=FALSE)

# Plot
p <- ggplot(data=rank_abund, aes(x=Rank, y=Abundance, colour=is.neg)) +
  geom_point() +
  scale_y_continuous(trans="log10") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=c(`TRUE`="grey", `FALSE`="darkgreen"), name = "Is control/blank") +
  ylab("log-transformed relative abundance") + xlab("Abundance rank in sample")

ggsave(file=file.path(subdir, "phy_sp_rank_abundance_plot.png"), p, width=8, height=6)

#### Taxonomic richness and number of reads, with histogram margines
p = ggplot(phy_sp@sam_data, aes(x=unmapped_count, y=taxa_raw, color = is.neg)) +
  geom_point(size=3) +
  labs(x="Number of reads", y="Number of OTUs") +
  scale_color_manual(values=c(`TRUE`="grey", `FALSE`="darkgreen"), name = "Is control/blank") +
  scale_x_continuous(trans="log10") +
  theme(legend.position = "bottom")

p <- ggMarginal(p, type="histogram", size=2, groupFill=TRUE)

ggsave(file.path(subdir, "phy_sp_taxa_and_read_dist.png"), p, width=8, height=8)

#### NMDS plot
ord <- ordinate(phy_sp_norm, "PCoA", "euclidean")

# Samples
p <- plot_ordination(phy_sp_norm, ord, color="Order_grouped", , shape="is.neg", title="PCoA plot of samples") +
  theme(legend.position = "bottom") +
  scale_shape_manual(values=c(`TRUE`=0, `FALSE`=16), name = "Is control/blank") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left")

ggsave(file.path(subdir, "phy_sp_sample_PCoA.png"), p, width=8, height=6)

# Heatmap

p <- plot_taxa_heatmap(phy_sp, subset.top=100, transformation="clr",
                      VariableA=c("Species", "is.neg"),
                      annotation_colors = list("Species" = species_palette, "is.neg" = c(`TRUE`="black", `FALSE`="white")))

ggsave(file.path(subdir, "phy_sp_heatmap.png"), p, width=8, height=6)