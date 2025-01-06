##### CODIV #####

#### Running MAG-host codiversification tests

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(tibble)
library(paco)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script

source(file.path("..", "scripts", "plot_setup.R"))

setwd(subdir)

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

############################
#### TIDY UP HOST TREES ####
############################

# Get consensus tree and fix tip labels
host_consensus <- consensus.edges(host_trees, method="least.squares")

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Equus_quagga", replacement="Equus_burchellii") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_byronia") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

# Drop tips from domesticated species
domesticated <- c("Ovis aries", "Sus scrofa domesticus", "Equus caballus")
host_consensus <- drop.tip(host_consensus, domesticated)

# Do the same with the metadata and mag trees
domesticate_mags <- append(bac_meta %>% filter(Species %in% domesticated) %>% pull(label),
                           ar_meta %>% filter(Species %in% domesticated) %>% pull(label))

bac_meta <- bac_meta %>% filter(!label %in% domesticate_mags)
ar_meta <- ar_meta %>% filter(!label %in% domesticate_mags)

bac_tree <- drop.tip(bac_tree, domesticate_mags)
ar_tree <- drop.tip(ar_tree, domesticate_mags)

#################################
#### CODIVERSIFICATION TESTS ####
#################################

# Select MAG clades to test
select_mag_clades <- function(tree, min_tips, max_depth) {
  cat(paste("Finding clades with at least", min_tips, "tips and a maximum depth of", max_depth, "\n"))
  # Empty list to store nodes
  nodes_to_test <- c()
  # Get all nodes
  nodes <- tree$node.label
  for (node in tree$node.label) {
    subtree <- extract.clade(tree, node)
    # Check if the node has more than min_tips tips
    if (length(subtree$tip.label) >= min_tips & max(node.depth(subtree)) <= max_depth) {
      # Add to list
      nodes_to_test <- append(nodes_to_test, node)
    }
  }
  return(nodes_to_test)
}

# Get PACo input: the subset of the MAG tree, ther corresponding host tree and the links
get_paco_input <- function(node, mag_tree, mag_meta, host_tree) {
  # Get subtree
  mag_subtree <- extract.clade(mag_tree, node)
  # Get bins and their host species based on keywords
  links <- mag_meta %>% select(label, Species) %>% filter(label %in% mag_subtree$tip.label)
  # Subset_host_tree
  hosts.to.drop <- host_tree$tip.label[!(host_tree$tip.label %in% links$Species)]
  host_tree_filtered <- drop.tip(host_tree, hosts.to.drop)
  return(list(mag_subtree, host_tree_filtered, links))
}

# Run PACo
run_paco <- function(mag_tree, host_tree, links) {
   # Get distance matrices
  host_dist <- cophenetic(host_tree)
  mag_dist <- cophenetic(mag_tree)
  # Get presence-absence matrix
  mat <- links %>%
      mutate(presence = 1) %>%
      select(Species, label, presence)
  mat <- mat %>% pivot_wider(values_from = "presence", names_from = "label") %>% as.data.frame
  mat[is.na(mat)] <- 0
  mat <- mat %>% column_to_rownames("Species") %>% as.matrix
  # Prepare paco data
  paco_dat <- prepare_paco_data(host_dist, mag_dist, mat)
  paco_coord <- add_pcoord(paco_dat)
  # Run paco
  res <- PACo(paco_coord, nperm = 1000, seed = NA, method = "r0", symmetric = FALSE, proc.warnings = TRUE, shuffled = FALSE)
  res_table <- data.frame(correction = res$correction,
                          p.value = res$gof$p,
                          ss = res$gof$ss)
  return(res_table)
}

######################
#### RUN ANALYSIS ####
######################

analysis <- function(mag_tree, mag_meta, host_tree, nodes_to_test) {
  paco_results <- data.frame()
  skipped <- data.frame()
  
  n=1
  for (node in nodes_to_test) {
    cat(paste0(n, ":", length(nodes_to_test), " - node ", node, "\n"))
    # Get PACo input
    paco_input <- get_paco_input(node, mag_tree, mag_meta, host_tree)
    mag_tree <- paco_input[[1]]
    host_tree <- paco_input[[2]]
    links <- paco_input[[3]]
    # Get some info about the clade
    info <- data.frame(node = node,
              mags_tips = length(mag_tree$tip.label),
              host_tips = length(host_tree$tip.label))
    
    if (info$host_tips <= 2) {
      cat("Skipping", node, "with less or equal to 2 host tips\n")
      skipped <- rbind(skipped, info)
      n <- n + 1
      next
    }
    # Run PACo
    res <- run_paco(mag_tree, host_tree, links)
    res <- cbind(info, res)
    # Append to table
    paco_results <- rbind(paco_results, res)
    n <- n + 1
  }
  # Account for multiple testing
  paco_results$p.adjust <- p.adjust(paco_results$p.value, method="BH")
  return(paco_results)
}

### For Bacteria ####

nodes_to_test <- select_mag_clades(bac_tree, min_tips=5, max_depth=20)

bac_paco_results <- analysis(bac_tree, bac_meta, host_consensus, nodes_to_test)

write.csv(bac_paco_results, file=file.path(subdir, "bac_paco_results.csv"), row.names=FALSE)
write.csv(skipped, file=file.path(subdir, "bac_skipped.csv"), row.names=FALSE)

### For Archaea ####

nodes_to_test <- select_mag_clades(ar_tree, min_tips=5, max_depth=20)

ar_paco_results <- analysis(ar_tree, ar_meta, host_consensus, nodes_to_test)

write.csv(ar_paco_results, file=file.path(subdir, "ar_paco_results.csv"), row.names=FALSE)
write.csv(skipped, file=file.path(subdir, "ar_skipped.csv"), row.names=FALSE)

######################
#### PLOT RESULTS ####
######################

# Plot cophylogenies
cophyloplot <- function(host_tree, mag_tree, links, host_colours) {
  links_rename <- links %>% rename("phy1"="label", "phy2"="Species")
  # Add colours
  links_rename$colour <- host_colours[match(links_rename$phy2, names(host_colours))]
  links_rename <- links_rename %>% as.matrix
  # Plot comparison
  coph <- cophylo(tr1=mag_tree, tr2=host_tree, assoc=links_rename)
  # Plot
  return(coph)
}

create_plots <- function(mag_tree, mag_meta, host_tree, paco_results, outdir) {
  n=1
  for (i in 1:nrow(paco_results)) {
    node <- paco_results$node[i]
    cat(paste0(n, ":", nrow(paco_results), " - node ", node, "\n"))
    # Get PACo input again: same files used in cophylo
    paco_input <- get_paco_input(node, mag_tree, mag_meta, host_tree)
    mag_tree <- paco_input[[1]]
    host_tree <- paco_input[[2]]
    links <- paco_input[[3]]
    # Plot cophyloplo
    coph <- cophyloplot(host_tree=host_tree, mag_tree=mag_tree, links=links, host_colours=species_palette)
    # Save plot
    png(file.path(subdir, outdir, paste(node, "cophyloplot.png", sep="_")), width = 1200, height = 600)
    plot(coph, link.type="curved",link.lwd=4, link.lty="solid", link.col=coph$assoc[,"colour"])
    text(x = 0.5, y = 0.95, labels = paste("adjusted p-value:", format(paco_results$p.adjust[i], digits = 3)), pos = 2, cex = 1.5, col = "black")
    dev.off()
    n <- n + 1
  }
}


# Save plots
dir.create(file.path(subdir, "bac_cophylo_plots"), showWarnings = FALSE)

create_plots(bac_tree, bac_meta, host_consensus, bac_paco_results, "bac_cophylo_plots")

# Save plots
dir.create(file.path(subdir, "ar_cophylo_plots"), showWarnings = FALSE)

create_plots(ar_tree, ar_meta, host_consensus, ar_paco_results, "ar_cophylo_plots")
