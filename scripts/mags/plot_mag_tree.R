##### Plot MAG tree #####

#### Plot MAG trees for Bacteria and Archaea and annotate using relevant metadata

#### LOAD PACKAGES ####
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# Bin metadata
bac_meta <- read.table(file.path(subdir, "bac_meta.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(subdir, "ar_meta.tsv"), sep="\t", header=TRUE)

# MAG trees
bac_tree <- read.tree(file = file.path(subdir, "bac_tree.tree"))
ar_tree <- read.tree(file = file.path(subdir, "ar_tree.tree"))

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

####################
#### PLOT TREES ####
####################

## Create phylum palettes
## Bacteria
bac_phyla <- bac_meta$phylum %>% unique
bac_phyla_df <- data.frame(phylum=bac_phyla,
                            # Colour all subgroups of Bacillota the same
                            grouped = ifelse(grepl("Bacillota_", bac_phyla), "Bacillota", bac_phyla))

# Create a custom palette using RColorBrewer
mag_palette_b <- left_join(bac_phyla_df,
                            data.frame(colour = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(bac_phyla_df$grouped))),
                             grouped = unique(bac_phyla_df$grouped)))

mag_palette_b <- setNames(mag_palette_b$colour, mag_palette_b$phylum)             

## Archeae
ar_phyla <- ar_meta$phylum %>% unique

mag_palette_a <- brewer.pal(length(ar_phyla), "Accent")
names(mag_palette_a) <- ar_phyla

#### Bacteria tree ####
# Colour by order and habitat
bac_p <- ggtree(bac_tree, layout = "circular", aes(color=phylum)) %<+% bac_meta +
  scale_colour_manual(values = mag_palette_b, name = "Phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=1.5, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=1, stroke=0.5, 
                aes(fill=host_order, color=habitat.general)) +
  scale_fill_manual(values = order_palette, name = "Order", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Habitat", na.value = "black") +
  scale_x_continuous(expand = c(0.1, 0.1)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(-10, -10, -10, -10), "cm"), # Remove margins
  legend.position=c(0.2, 0.7),
  legend.text = element_text(size=20),
  legend.title = element_text(size=20)) +
  guides(fill = guide_legend(override.aes = list(size = 5)), 
  color = guide_legend(override.aes = list(size = 5)))

# Highlight the phyla
bac_node_colours <- data.frame(node=integer(), phylum=character())
for (val in unique(bac_meta$phylum)) {
    if (is.na(val)) {
        next
    }
    tips <- bac_meta$label[bac_meta$is.tip & bac_meta$phylum == val]
    node <- ifelse(length(tips) > 1, getMRCA(bac_tree, tips), bac_meta$node[bac_meta$label == tips])
    bac_node_colours <- rbind(bac_node_colours, data.frame(node=node, phylum=val))
}
bac_node_colours$colour <- mag_palette_b[bac_node_colours$phylum]

for (i in 1:nrow(bac_node_colours)) {
    bac_p <- bac_p + geom_hilight(node=bac_node_colours$node[i], label=bac_node_colours$phylum[i], geom="label",
                                     fill=bac_node_colours$colour[i], alpha=0.2)
}

#### Archaea tree ####
# Colour by order and habitat
ar_p <- ggtree(ar_tree, layout = "circular", aes(color=phylum)) %<+% ar_meta +
  scale_colour_manual(values = mag_palette_a, name = "Phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=habitat.general)) +
  scale_fill_manual(values = order_palette, name = "Order", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Habitat", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(legend.position=c(0.1, 0.85),
  legend.text = element_text(size=10),
  legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
  color = guide_legend(override.aes = list(size = 2.5)))


# Highlight the phyla
ar_node_colours <- data.frame(node=integer(), phylum=character())
for (val in unique(ar_meta$phylum)) {
    if (is.na(val)) {
        next
    }
    tips <- ar_meta$label[ar_meta$is.tip & ar_meta$phylum == val]
    node <- ifelse(length(tips) > 1, getMRCA(ar_tree, tips), ar_meta$node[ar_meta$label == tips])
    ar_node_colours <- rbind(ar_node_colours, data.frame(node=node, phylum=val))
}
ar_node_colours$colour <- mag_palette_a[ar_node_colours$phylum]

for (i in 1:nrow(ar_node_colours)) {
    ar_p <- ar_p + geom_hilight(node=ar_node_colours$node[i], label=ar_node_colours$phylum[i], geom="label",
                                     fill=ar_node_colours$colour[i], alpha=0.2)
}

#####################
#### SAVE OUTPUT ####
#####################

ggsave(bac_p, file=file.path(subdir, "bac_genome_tree.png"), width = 20, height = 20)
ggsave(ar_p, file=file.path(subdir, "ar_genome_tree.png"), width = 12, height = 12)
