##### TEST HYPOTHESIS OF MICROBIOME EVOLUTION #####
## This script tests our hypotheses regarding the evolution of the microbiome in the context of the host phylogeny and ecology
## using statistical tests and producing relevant plots.

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(ape)
library(phytools)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
#library(microbiomeutilities)
#library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with output from the community analysis
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "hypothesis_testing")) # subdirectory for the output of this script
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

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

# Microbial taxa known uses
taxa_uses <- read.csv(file.path(outdir, "use_relations.csv"), header=TRUE)

########################
#### PHYLOSYMBIOSIS ####
########################

# Keep only wild animals
phy_wild <- phy_sp_f %>% subset_samples(!(Species %in% c("Ovis aries", "Equus caballus", "Sus scrofa domesticus")))
phy_wild <- phy_wild %>% subset_taxa(taxa_sums(phy_wild) > 0)

phy_wild_clr <- transform(phy_wild, "clr")

## Do host distances correlate with microbiome distances?

# Get host distances
# Get consensus tree and fix tip labels
host_consensus <- consensus.edges(host_trees, method="least.squares")

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Equus_quagga", replacement="Equus_burchellii") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_byronia") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

# Species distances
host_dist_melt <- cophenetic(host_consensus) %>% as.data.frame() %>% rownames_to_column("Item1") %>%
                melt(idvar = "Item1") %>% rename(Item2 = variable, host_distance = value)

# Samples to species
samples_to_species <- phy_wild@sam_data %>% data.frame %>% select(Species) %>% rownames_to_column("Sample")

# Add to dist
samples_dist_melt <- host_dist_melt %>%
                    right_join(samples_to_species, by=c("Item1"="Species"), relationship = "many-to-many") %>%
                    right_join(samples_to_species, by=c("Item2"="Species"), relationship = "many-to-many", suffix = c("1","2")) %>%
                    select(Sample1, Sample2, host_distance)

samples_dist <- dcast(samples_dist_melt, Sample1 ~ Sample2, value = "host_distance") %>% column_to_rownames("Sample1") %>% as.dist

# Microbiome distances
mb_dist <- vegdist(t(phy_wild_clr@otu_table), method="euclidean")

## Mantel test
mantel(samples_dist, mb_dist, method="pearson", permutations=999)

## Plot

# Melt microbiome distances
mb_dist_melt <- mb_dist %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample1") %>%
                melt(idvar = "Sample1") %>% rename(Sample2 = variable, microbiome_distance = value)

# Merge with host distances
dist_df <- samples_dist_melt %>% full_join(mb_dist_melt, by=c("Sample1"="Sample1", "Sample2"="Sample2"))

# Add info on diet and habitat
eco_info <- phy_wild@sam_data %>% data.frame %>% select(calculated_species_main_diet, habitat.general) %>%
                rownames_to_column("Sample") %>% rename(host_diet = calculated_species_main_diet,
                                                        host_habitat = habitat.general)

# Add to dist df 
dist_df <- dist_df %>% left_join(eco_info, by=c("Sample1"="Sample")) %>%
            left_join(eco_info, by=c("Sample2"="Sample"), suffix = c("1","2"))  %>%
            mutate(same_diet = ifelse(host_diet1 == host_diet2, "Same", "Different"),
                   same_habitat = ifelse(host_habitat1 == host_habitat2, "Same", "Different")) %>%
            # Remove intraspecific comparisons
            filter(host_distance > 0)

p <- ggplot(dist_df, aes(x=host_distance, y=microbiome_distance, colour = same_diet, shape = same_habitat)) +
  geom_jitter(width = 5, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = c("Same" = "#021020", "Different" = "#558CD2")) +
  scale_shape_manual(values = c("Same" = 16, "Different" = 3)) +
  geom_smooth(method="lm", aes(group = same_diet)) +
  labs(x="Host distance", y="Microbiome distance") +
  theme(legend.position="left")

ggsave(file.path(subdir, "host_microbiome_distance_correlation.png"), p, width=8, height=8)

# The same but only within diets
dist_df_same_diet <- dist_df %>% filter(same_diet == "Same")

p <- ggplot(dist_df_same_diet, aes(x=host_distance, y=microbiome_distance, colour = host_diet1, shape = same_habitat)) +
  geom_jitter(width = 5, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = diet_palette) +
  scale_shape_manual(values = c("Same" = 16, "Different" = 3)) +
  geom_smooth(method="lm", aes(group = host_diet1)) +
  labs(x="Host distance", y="Microbiome distance") +
  theme(legend.position="left")

ggsave(file.path(subdir, "host_microbiome_distance_correlation_same_diet.png"), p, width=5, height=5)

###################################################
#### MICROBIAL PHENOTYPES ASSOCIATED WITH DIET ####
###################################################

taxa_uses_wide <- taxa_uses %>% mutate(OBT = gsub(OBT, pattern = " ", replacement = "_")) %>% pivot_wider(id_cols = taxon, names_from = OBT, values_from = occurences) %>%
          rename("fermentative_activity" = "fermentative") %>%
          # Turn all values more than 1 to 1 and all NAs to 0
          mutate(across(ends_with("activity"), ~ifelse(is.na(.), 0, ifelse(. > 0, 1, 0)))) %>% 
          mutate(across(ends_with("activity"), as.numeric))

# Combine phyloseq info with taxa uses
taxa_uses_abund <- psmelt(transform(phy_sp_f, "compositional")) %>% select(Sample, OTU, Abundance, Species, Common.name, Order, calculated_species_main_diet, nfe, cf, cp, ee) %>%
                      left_join(taxa_uses_wide, by=c("OTU"="taxon"))

## Does the abundance of taxa with unknown use vary by species?
known_uses_abund <- psmelt(transform(phy_sp_f, "compositional")) %>% select(Sample, OTU, Abundance, Order, Common.name, calculated_species_main_diet, nfe, cf, cp, ee) %>%
                      filter(OTU %in% taxa_uses$taxon) %>% group_by(Sample, Order, Common.name, calculated_species_main_diet) %>% summarise(known_abund = sum(Abundance)) %>%
                      # Shorten "Omnivore" and "Animalivore" for plotting
                      mutate(calculated_species_main_diet = recode(calculated_species_main_diet, "Omnivore" = "O.", "Animalivore" = "Anim."))

# Get average known-use abundance across samples
mean_abund <- mean(known_uses_abund$known_abund)

p <- ggplot(known_uses_abund, aes(y = known_abund, x = Common.name, fill = Order)) +
  facet_grid(cols=vars(calculated_species_main_diet), space = "free", scales="free") +
  geom_violin(aes(colour = Order)) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  geom_hline(yintercept = mean_abund, linetype = "dotted") +
  scale_fill_manual(values = order_palette) +
  scale_colour_manual(values = order_palette) +
  labs(y="Abundance of taxa with\nat least on reported use") +
  theme(legend.position="bottom", axis.text.x = element_text(hjust = 1))

ggsave(file.path(subdir, "known_function_abundance.png"), p, width=8, height=7)

# Because species differ to the available info about the taxa they contain,
# I will normalise that abundance by the total abundance of the taxa in the sample
taxa_uses_abund_norm <- taxa_uses_abund %>% group_by(Sample) %>% filter(OTU %in% taxa_uses$taxon) %>% mutate(Total_abundance = sum(Abundance)) %>%
                        mutate(Norm_abundance = Abundance / Total_abundance) %>% ungroup()

## Does the abundance and diversity of proteolytic microorganisms correlate with the amount of crude protein in the diet?
proteolytic_abund <- taxa_uses_abund_norm %>% filter(proteolytic_activity > 0) %>% group_by(Sample, Species, Order, cp) %>% summarise(proteolytic_abund = sum(Norm_abundance))

p <- ggplot(proteolytic_abund, aes(y = proteolytic_abund, x = cp, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = order_palette) +
  labs(x="Crude protein in diet", y="Rel. abundance of proteolytic microorganisms (normalised)")

ggsave(file.path(subdir, "proteolytic_abundance_cp.png"), p, width=8, height=8)

## Does the abundance and diversity of fermentative microorganisms correlate with the amount of carbohydrates in the diet?
fermentative_abund <- taxa_uses_abund_norm %>% filter(fermentative_activity > 0) %>% group_by(Sample, Species, Order, nfe) %>% summarise(fermentative_abund = sum(Norm_abundance))

p <- ggplot(fermentative_abund, aes(y = fermentative_abund, x = nfe, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = order_palette) +
  labs(x="Nitrogen free extract in diet", y="Abundance of fermentative microorganisms")

ggsave(file.path(subdir, "fermentative_abundance_cp.png"), p, width=8, height=8)

## Does the abundance and diversity of lipolytic microorganisms correlate with the amount of ether extract in the diet?
lipolytic_abund <- taxa_uses_abund_norm %>% filter(lipolytic_activity > 0) %>% group_by(Sample, Species, Order, ee) %>% summarise(lipolytic_abund = sum(Norm_abundance))

p <- ggplot(lipolytic_abund, aes(y = lipolytic_abund, x = ee, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = order_palette) +
  labs(x="Ether extract in diet", y="Abundance of lipolytic microorganisms")

ggsave(file.path(subdir, "lipolytic_abundance_cp.png"), p, width=8, height=8)

## Control: antimicrobial activity and diet
antimicrobial_abund <- taxa_uses_abund_norm %>% filter(antimicrobial_activity > 0) %>% group_by(Sample, Species, Order, cp, nfe, ee) %>% summarise(antimicrobial_abund = sum(Norm_abundance))

p <- ggplot(antimicrobial_abund, aes(y = antimicrobial_abund, x = cp, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = order_palette) +
  labs(x="Crude protein in diet", y="Abundance of antimicrobial microorganisms")

ggsave(file.path(subdir, "antimicrobial_abundance_cp.png"), p, width=8, height=8)

ntaxa_uses <-
    taxa_uses_abund_norm %>% group_by(Sample, Common.name, Order, calculated_species_main_diet) %>%
      summarise(p_antimicrobial = sum(antimicrobial_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_antimicrobial = sum(Abundance[antimicrobial_activity > 0]),
                p_fermentative = sum(fermentative_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_fermentative = sum(Abundance[fermentative_activity > 0]),
                p_proteolytic = sum(proteolytic_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_proteolytic = sum(Abundance[proteolytic_activity > 0]),
                p_lipolytic = sum(lipolytic_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_lipolytic = sum(Abundance[lipolytic_activity > 0])) %>% 
                pivot_longer(cols = starts_with("p_") | starts_with("a_"), names_to = "activity", values_to = "value") %>%
                separate(activity, into = c("measure", "activity"), sep = "_", remove = FALSE) %>%
                mutate(measure = ifelse(measure == "p", "Proportion of taxa", "Proportion of abundance"))

# Summarise by species
ntaxa_uses_species <- 
    ntaxa_uses %>% group_by(Common.name, activity, measure, calculated_species_main_diet, Order) %>%
      summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE)) %>%
      ungroup() %>%
      # Shorten "Omnivore" for plotting
      mutate(calculated_species_main_diet = recode(calculated_species_main_diet, "Omnivore" = "Omn.", "Animalivore" = "Anim."))

p <- ggplot(ntaxa_uses_species, aes(x = Common.name, y = mean, fill = activity)) +
        geom_bar(stat = "identity", position = "dodge")  +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), colour = "black", position = position_dodge(width = 0.9), width = 0.25) +
        scale_fill_manual(values = c("antimicrobial" = "gray50",
                                     "fermentative" = "#3A459C",
                                     "proteolytic" = "#B41F34",
                                     "lipolytic" = "#BD9120")) +
        facet_grid(cols = vars(calculated_species_main_diet),
                   rows = vars(measure),
                   scales = "free", space = "free") +
        theme(legend.position = "bottom", axis.title = element_blank(),
              axis.text.x = element_text(hjust = 1)) 

ggsave(file.path(subdir, "activity_proportions_per_species.png"), p, width=12, height=8)
