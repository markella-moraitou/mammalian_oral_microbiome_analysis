##### PLOT SAMPLE METADATA #####

#### Summarise information about the samples and diet ####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(ape)
library(ggtree)
library(phytools)
library(rphylopic)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "sample_summary_plot")) # subdirectory for the output of this script
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

phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))
metadata <- data.frame(phy_sp_f@sam_data)

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

# Phylopic data for plotting
phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"))

########################
#### SUMMARISE DIET ####
########################

meta <- metadata %>% filter(!is.neg)

#### Diet PCA ####
diet_data <- unique(meta[,c("Common.name", "Species", "Order_grouped", "cp", "ee", "cf", "ash", "nfe", "diet.general")]) %>%
              filter(!grepl("blank", Species) & !grepl("control", Species))
              
rownames(diet_data) <- NULL
diet_data <- diet_data %>% column_to_rownames("Common.name")

diet_var <- diet_data[,c("cp", "ee", "cf", "ash", "nfe")]
ord <- prcomp(diet_var)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings_matrix <- data.frame(Variables = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")])

labels <- ord$x %>% as.data.frame %>% mutate(label = case_when(PC1 > 10 ~ rownames(.),
                                                               PC2 == max(PC2) | PC2 == min(PC2) ~ rownames(.))) %>%
  pull(label)

#labels <- ord$x %>% as.data.frame %>% mutate(label = rownames(.)) %>% pull(label)

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=diet_data$diet.general), data = data.frame(ord$x)) +
  geom_point(size=3) +
  scale_colour_manual(values = diet_palette, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings_matrix, aes(x = 0, y = 0, xend = (PC1*18),
                                       yend = (PC2*18)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  annotate("text", x = (loadings_matrix$PC1*20), y = (loadings_matrix$PC2*20),
           label = loadings_matrix$Variables) +
  xlim(c(-35, 80)) # Adjust position of text a bit

## Add principal components to data
meta$PC1 <- ord$x[match(meta$Common.name, rownames(ord$x)), 1]
meta$PC2 <- ord$x[match(meta$Common.name, rownames(ord$x)), 2]

ggsave(filename  =  file.path(subdir, "dietary_PCA.png"), pca, width  =  4, height = 4)

#### Diet summary ####
diet_long <- diet_data %>% rownames_to_column("Common.name") %>% pivot_longer(c(cp, ee, cf, ash, nfe), values_to = "proportion", names_to = "nutrient") %>%
  mutate(nutrient = factor(nutrient, levels = c("ash", "ee", "cp", "nfe", "cf"))) %>%
  mutate(diet.general = recode(diet.general, "Omnivore" = "Om.", "Animalivore" = "Anim.")) %>%
  arrange(nutrient)

nutrient_palette <- c(ash = "grey",
                      cf = "#3BA01B",
                      nfe = "#3A459C",
                      cp = "#B41F34",
                      ee = "#BD9120")

diet_barplot <- ggplot(diet_long, aes(x = Common.name, y = proportion, fill = nutrient, group = diet.general)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       legend.position = "right") +
  facet_grid(cols = vars(diet.general), scales = "free", space = "free") +
  scale_fill_manual(values = nutrient_palette,
                    labels = c("ash" = "ash\n(inorganic)",
                               "cf" = "CF\n(crude fiber)",
                               "nfe" = "NFE\n(carbohydrates)",
                               "cp" = "CP\n(crude protein)",
                               "ee" = "EE\n(fat)"),
                    name = "")

ggsave(filename  =  file.path(subdir, "diet_barplot.png"), diet_barplot, width  =  8, height = 4)

#### Order and diet ####
summ_dat <- meta %>% group_by(Order, diet.general) %>% 
  summarise(n_samples = n_distinct(new_name), n_species = n_distinct(Species))

write.csv(summ_dat, file = file.path(subdir, "data_summary.csv"), quote  =  FALSE, row.names  =  FALSE)

sums <- summ_dat %>% ungroup %>% summarise(total_species = sum(n_species), total_samples = sum(n_samples))

p_summary <-
  ggplot(aes(x = diet.general, y = Order, size = n_species, fill = n_samples, colour = n_samples), data = summ_dat) +
  geom_point(shape = 21) +
  scale_fill_continuous(low = "#CDB139", high = "#C73842", name = "Sample number", trans = "log", breaks = c(5, 25, 125)) +
  scale_colour_continuous(low = "#CDB139", high = "#C73842", name = "Sample number", trans = "log", breaks = c(5, 25, 125)) +
  scale_size(range  =  c(5, 18), breaks = c(min(summ_dat$n_species), mean(unique(summ_dat$n_species)), max(summ_dat$n_species)),
             name = "Species number") +
  new_scale(new_aes = "size") +
  geom_text(aes(label = paste0(n_species, "(", n_samples, ")"), x = diet.general, y = Order, size = n_species), fontface = "bold", colour = "black") +
  scale_size_continuous(range = c(2, 3), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(caption = "Number of species (number of samples)")

ggsave(filename  =  file.path(subdir, "dataset_summary_plot.png"), p_summary, width  =  4, height = 4)

########################
#### HOST PHYLOGENY ####
########################

# Get consensus tree and fix tip labels
host_consensus <- consensus.edges(host_trees, method="least.squares")

host_consensus$edge.length <- NULL

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Equus_quagga", replacement="Equus_burchellii") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_byronia") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

# Drop tips not in dataset
host_consensus <- drop.tip(host_consensus, setdiff(host_consensus$tip.label, meta$Species))

# Get traits for tree tips
host_traits <- data.frame(tip = host_consensus$tip.label,
                          node=nodeid(host_consensus, host_consensus$tip.label),
                          Common.name = meta$Common.name[match(host_consensus$tip.label, meta$Species)],
                          Order = meta$Order[match(host_consensus$tip.label, meta$Species)],
                          diet.general = meta$diet.general[match(host_consensus$tip.label, meta$Species)],
                          habitat.general = meta$habitat.general[match(host_consensus$tip.label, meta$Species)],
                          phylopic = phylopics$uid[match(host_consensus$tip.label, phylopics$Species)])

# Create a function that will assign the traits of the tips to all parent nodes until the MRCA
tips_to_nodes <- function(tree, trait_table, trait) {
  # Create a table to store the info
  node_table <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(node_table) <- c("node", sym(trait))
  # For each unique trait value
  for (value in unique(trait_table[, trait])) {
    # Get tips
    tips <- trait_table %>% filter(!!sym(trait) == value) %>% pull(tip)
    # If there is only one tip, the trait doesn't get assigned to any parent nodes
    if (length(tips)==1) {
      node_table <- node_table %>% rbind(data.frame(node = tip_index <- which(tree$tip.label == tips),
                                                    trait_value = value) %>% 
                                           rename(!!trait := trait_value))
    } else {
      # Get their mrca
      anc <- getMRCA(tree, tips)
      # Get descendant nodes of the mrca
      desc <- getDescendants(tree, anc)
      # Append the node numbers and corresponding trait values to node_table
      node_table <- node_table %>% rbind(data.frame(node = desc,
                                                    trait_value = rep(value, length(desc))) %>% 
                                           rename(!!trait := trait_value))
    }
  }
  return(node_table)
}

# Add order trait to nodes
node_traits <- tips_to_nodes(host_consensus, host_traits, "Order") %>% left_join(host_traits) %>% arrange(node)

# Add traits to tree
host_consensus_traits <- full_join(host_consensus, node_traits, by = "node")

# Plot tree
p <- 
  ggtree(host_consensus_traits, aes(colour=Order), layout = "circular", open.angle = 180, size=1) + 
  geom_tiplab(size=4, nudge_x = 4.5, aes(label = Common.name)) + 
  geom_phylopic(aes(uuid=phylopic, fill=Order), width = 0.7, position=position_nudge(x=3)) +
  # Colour branches by order
  scale_color_manual(values = order_palette, name = "Taxonomic order", na.value = "black") +
  scale_fill_manual(values = order_palette, name = "Taxonomic order") +
  # Add points indicating diet and habitat
  new_scale_colour() +
  new_scale_fill() +
  geom_tippoint(shape=21, size=4, stroke=2, 
                aes(fill=diet.general, colour=habitat.general), position = position_nudge(x = 1)) +
  scale_fill_manual(values = diet_palette, name = "Diet", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Habitat", na.value = "black") +
  # Adjust theme
  theme_tree() +
  theme(plot.margin = margin(t=0, r=3, b=0, l=3, "cm"),
        legend.position = "none")

ggsave(filename =  file.path(subdir, "host_phylogeny.png"), p, width  =  8, height = 8)
