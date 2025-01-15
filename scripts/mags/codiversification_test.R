##### CODIVERSIFICATION #####

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

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

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
domesticate_mags <- append(bac_meta %>% filter(host_species %in% domesticated) %>% pull(label),
                           ar_meta %>% filter(host_species %in% domesticated) %>% pull(label))

bac_meta <- bac_meta %>% filter(!label %in% domesticate_mags)
ar_meta <- ar_meta %>% filter(!label %in% domesticate_mags)

bac_tree <- drop.tip(bac_tree, domesticate_mags)
ar_tree <- drop.tip(ar_tree, domesticate_mags)

#################################
#### CODIVERSIFICATION TESTS ####
#################################

# Select MAG clades to test
# Outputs lists, startimg with the most recent clade with at least 5 tips
# Then go back until the maximum depth is reached
select_mag_clades <- function(tree, min_tips, max_depth) {
  cat(paste("Finding clades with at least", min_tips, "tips and a maximum depth of", max_depth, "\n"))
  # Start by identifying clades with at least 5 tips
  # but exclude parent nodes if the child node is already in the list
  initial_nodes <- c()
  cat("Finding initial nodes\n")
  for (node in Ntip(tree) + (tree$Nnode:1)) {
    descendants <- getDescendants(tree, node)
    if (length(descendants) >= 5 & length(intersect(descendants, initial_nodes)) == 0) {
      initial_nodes <- append(node, initial_nodes)
    }
  }
  # Empty list to store nodes
  nodes_lists <- list()
  cat("Going back to maximum depth\n")
  # Traverse tree backwards from each initial node to the maximum depth, adding nodes to a lineage list
  for (node in initial_nodes) {
    initial_node <- node
    lineage <- c()
    depth <- node.depth(tree)[node]
    while (depth <= max_depth & depth < max(node.depth(tree))) {
      depth <- node.depth(tree)[node]
      # Add to list
      lineage <- append(lineage, node)
      # Get the parent node
      node <- getParent(tree, node)
      depth <- node.depth(tree)[node]
    }
    # Add to list using node labels
    node_label <- tree$node.label[initial_node - Ntip(tree)]
    nodes_lists[[node_label]] <- tree$node.label[lineage - Ntip(tree)]
  }
  return(nodes_lists)
}

# Get input for codiversification analysis: the subset of the MAG tree, ther corresponding host tree and the links
prepare_input <- function(node, mag_tree, mag_meta, host_tree) {
  # Get subtree
  mag_subtree <- extract.clade(mag_tree, node)
  # Get bins and their host species based on keywords
  links <- mag_meta %>% select(label, host_species) %>% filter(label %in% mag_subtree$tip.label)
  # Subset_host_tree
  hosts.to.drop <- host_tree$tip.label[!(host_tree$tip.label %in% links$host_species)]
  host_tree_filtered <- drop.tip(host_tree, hosts.to.drop)
  return(list(mag_subtree, host_tree_filtered, links))
}

# Get correlation statistics
dist_correlation <- function(mag_tree, host_tree, links, plot=FALSE) {
   # Get distance matrices
  host_dist <- cophenetic(host_tree)[links$host_species, links$host_species] %>% as.dist
  mag_dist <- cophenetic(mag_tree) %>% as.dist
  # Run Mantel test
  r <- cor.test(host_dist, mag_dist, method="pearson")
  if (plot) {
    plot <- ggplot(data.frame(x=host_dist, y=mag_dist), aes(x=x, y=y)) +
            geom_point() +
            labs(x="Host distance", y="MAG distance") +
            labs(title = mag_tree$node.label[1], subtitle = paste("r:", format(r$statistic, digits=3), "p:", format(r$p.value, digits=3)), size=5)
    ggsave(plot=plot, filename=file.path(intermediatedir, paste0(mag_tree$node.label[1], "_cor.png")), width=4, height=4)
  }
  return(list(r=r$statistic, p.value=r$p.value))
}

# Permute links
permute_tree_tips <- function(tree) {
  tree_perm <- tree
  tree_perm$tip.label <- sample(tree$tip.label)
  return(tree_perm)
}

## TO DO: Use p-value instead of r to select the best node
scan_lineage <- function(lineage, mag_tree, mag_meta, host_tree) {
  min_p <- 1
  best_node <- NULL
  # For each node in the lineage list
  for (node in lineage) {
    # Get input
    input <- prepare_input(node, mag_tree, mag_meta, host_tree)
    mag_tree_sub <- input[[1]]
    host_tree_sub <- input[[2]]
    links <- input[[3]]
    if (length(host_tree_sub$tip.label) <= 2) {
      cat("Skipping", node, "with less or equal to 2 host tips\n")
      next
    }
    # Calculate the correlation between the host and MAG distance matrices and permute to get p-value
    res <- permutation_test(mag_tree_sub, host_tree_sub, links, nperm=500)
    r <- res$r
    p <- res$p.value
    cat(node, "- r:", r, " - p:", p, ": #tips:", length(mag_tree_sub$tip.label), "\n")
    # Keep the node with the highest correlation along the lineage
    if (p < min_p) {
      min_p <- r
      best_node <- node
    }
  }
  return(best_node)
}

permutation_test <- function(mag_tree, host_tree, links, nperm) {
    # Calculate the correlation between the host and MAG distance matrices
    r <- dist_correlation(mag_tree, host_tree, links, plot=TRUE)$r[[1]]
    # Permute host species and calculate correlation
    permuted_r <- replicate(nperm, dist_correlation(mag_tree, permute_tree_tips(host_tree), links)$r[[1]]) %>% 
                  unlist
    # Calculate p-value
    p_value <- sum(permuted_r >= r) / length(permuted_r)
    # Plot
    plot <- ggplot(data.frame(r=permuted_r), aes(x=r)) +
          geom_histogram() +
          geom_vline(xintercept=r, color="red") +
          labs(title = mag_tree$node.label[1], subtitle=paste("r:", format(r, digits=3), "p:", format(p_value, digits=3)), size=5)
    ggsave(plot=plot, filename=file.path(intermediatedir, paste0(mag_tree$node.label[1], "_perm.png")), width=4, height=4)
    
    return(list(r=r, p.value=p_value))
}

######################
#### RUN ANALYSIS ####
######################

analysis <- function(mag_tree, mag_meta, host_tree, nodes_to_test, intermediatedir) {
  # Create directory for intermediate results
  intermediatedir <<- intermediatedir
  dir.create(file.path(subdir, intermediatedir), showWarnings = FALSE)
  n <- 1
  lineages <- length(nodes_to_test)
  best_nodes <- c()
  for (lineage_name in names(nodes_to_test)) {
    cat(paste0("Lineage ", n, ":", lineages, " - ", lineage_name, "\n"))
    lineage <- nodes_to_test[[lineage_name]]
    # Scan lineage and get the node with the highest correlation
    best_node <- scan_lineage(lineage, mag_tree, mag_meta, host_tree)
    best_nodes <- append(best_nodes, best_node)
    if (is.null(best_node)) {
      next
    }
    n <- n + 1
  }
  best_nodes <- unique(best_nodes)
  cat("Will test", length(best_nodes), "nodes\n")
  # Then test only best nodes
  n <- 1
  total <- length(best_nodes)
  test_results <- data.frame()
  for (node in best_nodes) {
    cat(paste0("Node ", n, ":", total, " - ", node, "\n"))
    # Prepare analysis input
    input <- prepare_input(node, mag_tree, mag_meta, host_tree)
    mag_tree_sub <- input[[1]]
    host_tree_sub <- input[[2]]
    links <- input[[3]]
    # Get some info about the clade
    info <- data.frame(node = node,
              mags_tips = length(mag_tree_sub$tip.label),
              host_tips = length(host_tree_sub$tip.label))
    # Run permutation test
    perm = permutation_test(mag_tree_sub, host_tree_sub, links, nperm=500)
    r <- perm$r
    p_value <- perm$p.value
    # Combine with node info
    res <- data.frame(r=r, p.value=p_value)
    res <- cbind(info, res)
    # Append to table
    test_results <- rbind(test_results, data.frame(res))
    n <- n + 1
  }
  # Account for multiple testing
  test_results$p.adjust <- p.adjust(test_results$p.value, method="BH")
  return(test_results)
}

### For Bacteria ####

nodes_to_test <- select_mag_clades(bac_tree, min_tips=5, max_depth=20)

bac_cod_results <- analysis(bac_tree, bac_meta, host_consensus, nodes_to_test, "bac_cod_intermediate")

write.csv(bac_cod_results, file=file.path(subdir, "bac_cod_results.csv"), row.names=FALSE)
write.csv(skipped, file=file.path(subdir, "bac_skipped.csv"), row.names=FALSE)

### For Archaea ####

nodes_to_test <- select_mag_clades(ar_tree, min_tips=5, max_depth=20)

ar_cod_results <- analysis(ar_tree, ar_meta, host_consensus, nodes_to_test, "ar_cod_intermediate")

write.csv(ar_cod_results, file=file.path(subdir, "ar_cod_results.csv"), row.names=FALSE)
write.csv(skipped, file=file.path(subdir, "ar_skipped.csv"), row.names=FALSE)

######################
#### PLOT RESULTS ####
######################

# Plot cophylogenies
cophyloplot <- function(host_tree, mag_tree, links, host_colours) {
  links_rename <- links %>% rename("phy1"="label", "phy2"="host_species")
  # Add colours
  links_rename$colour <- host_colours[match(links_rename$phy2, names(host_colours))]
  links_rename <- links_rename %>% as.matrix
  # Plot comparison
  coph <- cophylo(tr1=mag_tree, tr2=host_tree, assoc=links_rename)
  # Plot
  return(coph)
}

create_plots <- function(mag_tree, mag_meta, host_tree, cod_results, outdir) {
  n=1
  for (i in 1:nrow(cod_results)) {
    node <- cod_results$node[i]
    cat(paste0(n, ":", nrow(cod_results), " - node ", node, "\n"))
    # Get PACo input again: same files used in cophylo
    input <- prepare_input(node, mag_tree, mag_meta, host_tree)
    mag_tree_sub <- input[[1]]
    host_tree_sub <- input[[2]]
    links <- input[[3]]
    # Plot cophyloplo
    coph <- cophyloplot(host_tree=host_tree_sub, mag_tree=mag_tree_sub, links=links, host_colours=species_palette)
    # Save plot
    png(file.path(subdir, outdir, paste(node, "cophyloplot.png", sep="_")), width = 1200, height = 600)
    par(mar=c(5, 4, 4, 2) + 0.1)
    plot(coph, link.type="curved",link.lwd=4, link.lty="solid", link.col=coph$assoc[,"colour"])
    # Add title
    title(paste("Node:", node, cex.main=1.5))
    text(x = 0.5, y = 0.95, labels = paste("p-value:", format(cod_results$p.value[i], digits = 3)), pos = 2, cex = 1.5, col = "black")
    dev.off()
    n <- n + 1
  }
}

# Save plots
dir.create(file.path(subdir, "bac_cophylo_plots"), showWarnings = FALSE)

create_plots(bac_tree, bac_meta, host_consensus, bac_cod_results, "bac_cophylo_plots")

# Save plots
dir.create(file.path(subdir, "ar_cophylo_plots"), showWarnings = FALSE)

create_plots(ar_tree, ar_meta, host_consensus, ar_cod_results, "ar_cophylo_plots")
