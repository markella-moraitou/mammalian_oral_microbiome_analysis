##### Get MAG tree #####

#### Filter GTDB-Tk output to keep only bins, and combine with metadata

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(renv)
library(stringr)
library(tibble)
library(ape)
library(tidytree)
library(tibble)
library(treeio)
library(phytools)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(rentrez)

#### VARIABLES AND WORKING DIRECTORY ####

indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

# Create subdirectory for output
dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

# Get ENTREZ KEY for taxonomic information
options(ENTREZ_KEY = Sys.getenv("API_KEY"))

#######################
#####  LOAD INPUT #####
#######################

# GTDB-Tk de_novo output
bac_tree <- read.tree(file.path(indir, "gtdbtk.bac120.decorated_itol.tree"))
ar_tree <- read.tree(file.path(indir, "gtdbtk.ar53.decorated_itol.tree"))

# bin metadata
mag_meta <- read.csv(file.path(indir, "hq_bin_metadata.csv"), header=TRUE)
# Replace more than one instance of . with a single .
colnames(mag_meta) <- gsub("..", ".", colnames(mag_meta), fixed=TRUE) %>% gsub("...", ".", ., fixed=TRUE)

# Host metadata
host_taxonomy <- read.csv(file.path(indir, "host_taxonomy.csv"), header=TRUE)
host_habitat <- read.csv(file.path(indir, "species_habitats.csv"), header=TRUE)
host_diet <- read.csv(file.path(indir, "Lintulaakso_diet_filtered.csv"), header=TRUE)

#################
#### PROCESS ####
#################

### Combine host metadata ###
colnames(host_taxonomy) <- paste0("host_", colnames(host_taxonomy))

host_meta <- host_taxonomy[, colnames(host_taxonomy)!="host_species"] %>%
    left_join(host_habitat, by=c("host_museum.species"="Species")) %>%
    left_join(host_diet,  by=c("host_museum.species"="Species"))

# Get taxon label: the most specific taxonomic level that is not empty and is not a made up GTDB name
get_taxon_label <- function(row) {
    # For each rank, return the first non-empty value that does not contain numbers
    for (rank in c("species", "genus", "family", "order", "class", "phylum")) {
        if (!is.na(row[[rank]]) & !grepl("[0-9]", row[[rank]])) {
            return(row[[rank]])
        }
    }
}

## Extract taxonomic metadata
meta <- mag_meta %>% rename("host_species"="Species") %>%
    left_join(host_meta, by=c("host_species"="host_museum.species"), relationship = "many-to-one") %>%
    separate(classification, sep=";", into=c("domain", "phylum", "class", "order", "family", "genus", "species")) %>%
    mutate(across(domain:species, ~ case_when(str_detect(.x, "^[dpcofgs]__$") ~ NA, TRUE ~ .x))) %>%
    mutate(is.mag=TRUE) %>%
    # Get taxonomy label
    mutate(classification=apply(., 1, get_taxon_label)) %>%
    group_by(classification) %>% mutate(label=paste(gsub(" ", "_", classification), row_number(), sep="_")) %>%
    # Remove *__ prefix
    mutate(across(domain:species, ~ str_remove(.x, "^[dpcofgs]__"))) %>% ungroup()

## Get taxonomic ids
get_taxon_id <- function(taxon_name) {
  search_results <- entrez_search(db = "taxonomy", term = paste(taxon_name, "AND (Bacteria[Subtree] OR Archaea[Subtree]) "))
  if (length(search_results$ids) == 1) {
    return(search_results$ids[1])
  } else if (length(search_results$ids) == 0) {
    cat("No results for", taxon_name, "\n")
    return(NA)
  } else {
    cat("Multiple results for", taxon_name, "\n")
    return(NA)
  }
 }

# Get unique classifications
name_to_taxid = data.frame(name = unique(meta$classification))

# Get simplied classifications to facilitate searching the NCBI DB
# e.g. remove underscores after genus and species names, as well as the *__ prefix
name_to_taxid$simple_name <- name_to_taxid$name %>%
                            str_remove("^[dpcofgs]__") %>%
                            str_replace("_.+ ", " ") %>% # Remove underscore plus character(s) before a space
                            str_replace("_.+$", "") # Remove underscore plus character at the end of the string

# Get taxids                         
taxids <- sapply(name_to_taxid$simple_name, get_taxon_id)
name_to_taxid$taxids <- unlist(taxids)[match(name_to_taxid$simple_name, names(unlist(taxids)))]

# Save
write.csv(name_to_taxid, file = file.path(subdir, "name_to_taxid.csv"), row.names=FALSE, quote=FALSE)

# Add ids to metadata
meta$ncbi_taxid <- name_to_taxid$taxids[match(meta$classification, name_to_taxid$name)]

# Filter MAG metadata
bac_meta <- meta %>% filter(domain == "Bacteria")
ar_meta <- meta %>% filter(domain == "Archaea")

# Function to subset tree and add metadata
subset_tree <- function(tree, metadata) {
    # Get tree with only MAGs
    tips.to.drop <- tree$tip.label[!(tree$tip.label %in%  str_remove(meta$bin, ".gz"))]
    pruned_tree <- drop.tip(tree, tips.to.drop)
    # Rename labels
    pruned_tree$tip.label <- meta$label[match(pruned_tree$tip.label, str_remove(meta$bin, ".gz"))]
    #pruned_tree <- unroot(pruned_tree)
    # When the node label is missing, use the node number
    node_labs <- paste0("N", (Ntip(pruned_tree)+1):(Ntip(pruned_tree) + Nnode(pruned_tree)))
    pruned_tree$node.label <- ifelse(pruned_tree$node.label=="", node_labs, pruned_tree$node.label)
    return(pruned_tree)
}

# Bacteria tree
bac_tree_new <- subset_tree(bac_tree, bac_meta)

# Archaea tree
ar_tree_new <- subset_tree(ar_tree, ar_meta)

# Add some traits should be added to the entire clade that shares a common ancestor, not just the tips
add_node_metadata <- function(tree, metadata, traits) {
    # Add metadata to tree
    tree_tibble <- tree %>% as_tibble() %>% mutate(node=1:nrow(.))
    tree_tibble <- tree_tibble %>% left_join(metadata, by="label")
    tree_tibble$is.tip <- ifelse(tree_tibble$node < Ntip(tree), TRUE, FALSE)
    # Add metadata to nodes that share a common ancestor
    for (tr in traits) {
        node_traits <- data.frame(node=integer(), val=character())
        for (val in unique(tree_tibble[!is.na(tree_tibble[[tr]]), ][[tr]])) {
            # Get all tips with this value
            tips <- tree_tibble$label[!is.na(tree_tibble[[tr]]) & tree_tibble[[tr]] == val]
            # Get the node of the most recent common ancestor
            node <- ifelse(length(tips) > 1, getMRCA(tree, tips), tree_tibble$node[tree_tibble$label == tips])
            # Get the descendants of this node
            desc <- getDescendants(tree, node)
            node_traits <- rbind(node_traits, data.frame(node=desc, val=rep(val, length(desc))))
        }
        # Add to tree tibble
        tree_tibble[[tr]] <- node_traits$val[match(tree_tibble$node, node_traits$node)]
    }
    metadata_with_nodes <- tree_tibble %>% select(-parent, -node, -branch.length)
    return(metadata_with_nodes)
}

bac_meta_new <- add_node_metadata(bac_tree_new, bac_meta, c("phylum", "class", "order", "family", "genus"))
ar_meta_new <- add_node_metadata(ar_tree_new, ar_meta, c("phylum", "class", "order", "family", "genus"))

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
bac_p <- ggtree(bac_tree_new, layout = "circular", aes(color=phylum)) %<+% bac_meta_new +
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
for (val in unique(bac_meta_new$phylum)) {
    if (is.na(val)) {
        next
    }
    tips <- bac_meta_new$label[bac_meta_new$is.tip & bac_meta_new$phylum == val]
    node <- ifelse(length(tips) > 1, getMRCA(bac_tree_new, tips), bac_meta_new$node[bac_meta_new$label == tips])
    bac_node_colours <- rbind(bac_node_colours, data.frame(node=node, phylum=val))
}
bac_node_colours$colour <- mag_palette_b[bac_node_colours$phylum]

for (i in 1:nrow(bac_node_colours)) {
    bac_p <- bac_p + geom_hilight(node=bac_node_colours$node[i], label=bac_node_colours$phylum[i], geom="label",
                                     fill=bac_node_colours$colour[i], alpha=0.2)
}

#### Archaea tree ####
# Colour by order and habitat
ar_p <- ggtree(ar_tree_new, layout = "circular", aes(color=phylum)) %<+% ar_meta_new +
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
for (val in unique(ar_meta_new$phylum)) {
    if (is.na(val)) {
        next
    }
    tips <- ar_meta_new$label[ar_meta_new$is.tip & ar_meta_new$phylum == val]
    node <- ifelse(length(tips) > 1, getMRCA(ar_tree_new, tips), ar_meta_new$node[ar_meta_new$label == tips])
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

write.tree(bac_tree_new, file = file.path(subdir, "bac_tree.tree"))
write.tree(ar_tree_new, file = file.path(subdir, "ar_tree.tree"))

write.table(bac_meta_new, file = file.path(subdir, "bac_meta.tsv"), sep = "\t", row.names=FALSE, quote=FALSE)
write.table(ar_meta_new, file = file.path(subdir, "ar_meta.tsv"), sep = "\t", row.names=FALSE, quote=FALSE)

ggsave(bac_p, file=file.path(subdir, "bac_genome_tree.png"), width = 20, height = 20)
ggsave(ar_p, file=file.path(subdir, "ar_genome_tree.png"), width = 12, height = 12)
