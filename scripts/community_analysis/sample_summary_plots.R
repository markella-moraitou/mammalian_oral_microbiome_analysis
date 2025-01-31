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

########################
#### SUMMARISE DIET ####
########################

meta <- metadata %>% filter(!is.neg)

#### Diet PCA ####
diet_data <- unique(meta[,c("Common.name", "Species", "Order_grouped", "cp", "ee", "cf", "ash", "nfe", "calculated_species_main_diet")]) %>%
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
pca <- ggplot(aes(x = PC1, y = PC2, colour=diet_data$calculated_species_main_diet), data = data.frame(ord$x)) +
  geom_point(size=3) +
  scale_colour_manual(values = diet_palette, name = "Dietary category") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings_matrix, aes(x = 0, y = 0, xend = (PC1*18),
                                       yend = (PC2*18)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  theme(legend.position = "bottom") +
  annotate("text", x = (loadings_matrix$PC1*20), y = (loadings_matrix$PC2*20),
           label = loadings_matrix$Variables) +
  xlim(c(-35, 80)) + # Adjust position of text a bit
  labs(title = "PCA based on proximate analysis data",
       caption = "ash = inorganic, cf = crude fiber, cp = crude protein, ee = ether extract (proxy for fat),
       nfe = nitrogen-free extract (proxy for sugars and starches)")

## Add principal components to data
meta$PC1 <- ord$x[match(meta$Common.name, rownames(ord$x)), 1]
meta$PC2 <- ord$x[match(meta$Common.name, rownames(ord$x)), 2]

ggsave(filename  =  file.path(subdir, "dietary_PCA.png"), pca, width  =  8, height = 8)

#### Diet summary ####
diet_long <- diet_data %>% rownames_to_column("Common.name") %>% pivot_longer(c(cp, ee, cf, ash, nfe), values_to = "proportion", names_to = "nutrient") %>%
  mutate(nutrient = factor(nutrient, levels = c("ash", "ee", "cp", "nfe", "cf"))) %>%
  arrange(nutrient)

nutrient_palette <- c(ash = "grey",
                      cf = "#3BA01B",
                      nfe = "#3A459C",
                      cp = "#B41F34",
                      ee = "#BD9120")

diet_barplot <- ggplot(diet_long, aes(x = Common.name, y = proportion, fill = nutrient, group = calculated_species_main_diet)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") +
  facet_grid(cols = vars(calculated_species_main_diet), scales = "free", space = "free") +
  scale_fill_manual(values = nutrient_palette,
                    labels = c("ash" = "ash (inorganic)",
                               "cf" = "CF (crude fiber)",
                               "nfe" = "NFE (carbohydrates)",
                               "cp" = "CP (crude protein)",
                               "ee" = "EE (fat)"),
                    name = "Nutrient") +
  xlab("Species") +
  labs(title = "Composition of diet for study species",
       subtitle = "based on proximate analysis",
       tag = "A.")

ggsave(filename  =  file.path(subdir, "diet_barplot.png"), diet_barplot, width  =  12, height = 8)

#### Order and diet ####
summ_dat <- meta %>% group_by(Order, calculated_species_main_diet) %>% 
  summarise(n_samples = n_distinct(new_name), n_species = n_distinct(Species))

write.csv(summ_dat, file = file.path(subdir, "data_summary.csv"), quote  =  FALSE, row.names  =  FALSE)

sums <- summ_dat %>% ungroup %>% summarise(total_species = sum(n_species), total_samples = sum(n_samples))

# Add label colors
summ_dat <- summ_dat %>% mutate(lab_col = case_when(n_samples < 25 ~ "black", TRUE ~ "white"))

p_summary <-
  ggplot(aes(x = calculated_species_main_diet, y = Order, size = n_species, fill = n_samples), data = summ_dat) +
  geom_point(shape = 21, colour = "black") +
  xlab("Dietary category") + ylab("Taxonomic order") +
  labs(tag = "B.") +
  scale_fill_continuous(low = "lightblue", high = "darkblue", name = "Sample number", trans = "log", breaks = c(5, 25, 125)) +
  scale_size(range  =  c(12, 30), breaks = c(min(summ_dat$n_species), mean(unique(summ_dat$n_species)), max(summ_dat$n_species)),
             name = "Species number") +
  new_scale(new_aes = "size") +
  geom_text(aes(label = paste0(n_species, " (", n_samples, ")"), x = calculated_species_main_diet, y = Order, colour = lab_col, size = n_species), fontface = "bold") +
  scale_colour_manual(values = summ_dat$lab_col, guide = "none") +
  scale_size_continuous(range = c(3, 6), guide = "none")

ggsave(filename  =  file.path(subdir, "dataset_summary_plot.png"), p_summary, width  =  9, height = 8)

