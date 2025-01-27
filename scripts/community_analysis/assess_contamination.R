##### ASSESS CONTAMINANTION #####

#### Use various sources of information to assess if a taxon is a contaminant or not
#### and if a sample is contaminated or not

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(tidyr)
library(tibble)
library(scales)
library(decontam)
library(cuperdec)
library(cowplot)
library(microbiomeutilities)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # subdirectory for the output of this script

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp <- readRDS(file.path(subdir, "phyloseq_objects", "phy_sp.RDS"))
phy_sp_norm <- readRDS(file.path(subdir, "phyloseq_objects", "phy_sp_norm.RDS"))

# Habitat OTU relations
habitats_table <- read.csv(file.path(subdir, "habitat_relations.csv"))

# Set minimum sample content threshold
min_samp <- 10^5

# Set minimum relative abundance threshold
min_ab <- 10^-3

##############################
#### COLLECT INFO ON TAXA ####
##############################

# Use information regarding the abundance ratio between samples and controls,
# prevalence, mean abundance and habitat info to assess if a taxon is a contaminant or not

## Collect info that will help decide if a taxon is a contaminant or not
phy_sp@sam_data$sample_type <- factor(ifelse(!phy_sp@sam_data$is.neg, "sample",
                                             ifelse(phy_sp@sam_data$Species == "Environmental control", "swab", "blank")),
                                             levels = c("sample", "swab", "blank"))

#### RELATIVE ABUNDANCE IN SAMPLES VS CONTROLS ####
# Get relative abundances
phy_sp_r <- transform(phy_sp, "compositional")
phy_sp_m <- phy_sp_r %>% psmelt 

# Is a taxons abundance higher in samples or environmental controls?
# Where does a sample have its highest relative abundance? Sample, blank or control?
abundance_ratios <- 
            phy_sp_m %>%
            # Add pseudocount to avoid division by zero
            mutate(Abundance = Abundance + 1/max(unmapped_count, na.rm = TRUE)) %>%
            # Get mean and max abundance per OTU in samples, controls and blanks
            group_by(OTU) %>%
            mutate(sample_type = case_when(grepl("control", Species) ~ "control",
                                           grepl("blank", Species) ~ "blank",
                                           !is.neg ~ "sample")) %>%
            group_by(OTU) %>%
            summarise(mean_abundance_samples = mean(Abundance[sample_type == "sample"]),
                      mean_abundance_controls = mean(Abundance[sample_type == "control"]),
                      ) %>%
            ungroup %>%
            # Get ratios of mean and max abundances
            mutate(mean_ratio = mean_abundance_samples/mean_abundance_controls) %>% select(OTU, mean_ratio)

otu_order <- abundance_ratios %>% arrange(mean_ratio) %>% pull(OTU)

abundance_ratios$OTU <- factor(abundance_ratios$OTU, levels = otu_order)

abundance_ratios <- abundance_ratios %>% arrange(OTU)

# Draw threshold of OTUs double as abundant in samples
ythresh <- abundance_ratios %>% filter(mean_ratio > 2) %>% slice_min(mean_ratio, n = 1) %>% pull(OTU)

# Plot
p_a <- ggplot(abundance_ratios, aes(x = mean_ratio, y = OTU)) +
  geom_bar(stat = "identity") +
  theme(legend.position="top") +
  ylab("OTU") +
  scale_x_continuous(name = "\nmean rel. abund. ratio in\nsamples/negatives",
                    trans = "log10", breaks = c(1, 100)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = ythresh, linetype = "dashed") +
  theme(axis.text.y = element_blank(), axis.ticks.x = element_blank())

##### PREVALENCE AND AVERAGE RELATIVE ABUNDANCE ####
# Get prevalence of OTU in samples, controls and blanks
prevalence_df <-
  rownames_to_column(data.frame(prevalence = prevalence(subset_samples(phy_sp_r, !is.neg), detection = min_ab)), "OTU") %>% 
  mutate(OTU = factor(OTU, levels = otu_order))

p_p <- ggplot(prevalence_df, aes(x = prevalence, fill = prevalence, y = OTU)) +
  geom_point(shape = 21) +
  scale_colour_viridis_c(option = "magma") + xlab("\nprevalence") +
  geom_hline(yintercept = ythresh, linetype = "dashed") +
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), , axis.title.x.top = element_text())

abundance_df <-
  phy_sp_m %>% filter(!is.neg) %>% group_by(OTU) %>% summarise(mean_abundance = mean(Abundance)) %>% ungroup %>%
  mutate(OTU = factor(OTU, levels = otu_order), 
        # Add pseudocount to avoid division by zero
        mean_abundance = mean_abundance + min_ab/10)

# Plot
p_ma <- ggplot(abundance_df, aes(x = mean_abundance, fill = mean_abundance, y = OTU)) +
  geom_point(shape = 21) +
  scale_colour_viridis_c() +
  scale_x_log10() + xlab("mean\nrel. abundance") +
  # Add threshold line
  geom_vline(xintercept = min_ab, linetype = "dashed") +
  geom_hline(yintercept = ythresh, linetype = "dashed") +
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x.top = element_text())

#### TAXON HABITATS ####

# For each taxon, get the number of references to each habitat
habitats <- data.frame(taxon = otu_order) %>% left_join(habitats_table) %>%
            # Fill NAs with 0
            mutate(occurences = replace(occurences, is.na(occurences), 0)) %>%
            # Group OBTs into larger habitats
            mutate(habitat = case_when(OBT %in% c("dental plaque", "mouth") ~ "oral",
                                       OBT %in% c("marine water", "deep sea") ~ "marine",
                                       OBT %in% c("mammalian", "wild animal", "mammalian livestock") ~ "animal",
                                       TRUE ~ OBT)) %>%
            # Make sure that each taxon has an occurence value for each habitat, even if 0, instead of row missing
            tidyr::complete(taxon, habitat, fill = list(occurences = 0)) %>%
            group_by(taxon, habitat) %>%
            summarise(occurences = sum(occurences)) %>%
            mutate(habitat = gsub(" ", "_", habitat)) %>%
            # Only select some of the most informative terms
            filter(habitat %in% c("oral", "animal", "soil", "marine", "rumen", "gut"))

# Get occurences as a percentage of the total
habitats <- habitats %>% group_by(taxon) %>% mutate(perc_occurences = occurences/sum(occurences)) %>% ungroup

# Set habitat as a factor
habitats$habitat <- factor(habitats$habitat, levels = c("oral", "animal", "rumen", "gut", "marine", "soil", "laboratory_equipment"))

habitats$taxon <- factor(habitats$taxon, levels = otu_order)

#### Combine all tables ####
# Get a wider version of the habitats table
habitats_wide <- habitats %>% select(taxon, habitat, perc_occurences) %>% pivot_wider(names_from = "habitat", values_from = "perc_occurences")

# Combine sample/negative abundance ratio, prevalence, mean abundance and habitat info
assess_taxa <- abundance_ratios %>%
               left_join(prevalence_df) %>%
               left_join(abundance_df) %>%
               left_join(habitats_wide, by = c("OTU" = "taxon"))

# Save table
write.table(assess_taxa, file=file.path(subdir, "assess_taxa.csv"), sep=",", row.names=FALSE, quote=FALSE)

#### Plot ####
p_h <- habitats %>%
      # Replace 0 with NA for plotting
      mutate(perc_occurences = replace(perc_occurences, perc_occurences == 0, NA)) %>%
      ggplot(aes(x = habitat, colour = habitat, fill = habitat, y = taxon, size = perc_occurences)) +
      geom_point(shape = 21, alpha = 0.2) +
      scale_size_continuous(range = c(0.5, 5)) +
      geom_hline(yintercept = ythresh, linetype = "dashed") +
      scale_fill_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73")) +
      scale_color_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73")) +
      theme(, axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x.top = element_text())

p <- plot_grid(p_a + theme(legend.position="bottom"),
               p_p + theme(legend.position="none"),
               p_ma + theme(legend.position="none"),
               p_h + theme(legend.position="none"), nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 0.75, 0.75, 1))

ggsave(file=file.path(subdir, "assess_taxa.png"), p, width=8, height=16)

#################################
#### COLLECT INFO ON SAMPLES ####
#################################

# Use different sources of information to assess if a sample is contaminated or not

#### DECOM ####
# Get decom results
decom_tbl <- data.frame(phy_sp@sam_data) %>%
  select(new_name, Common.name, Species, Order_grouped, is.neg, starts_with("p_")) %>%
  # turn all NAs to 0 (for plotting)
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  # Calculate oral to soil+skin ratio
  mutate(oral_contam_ratio = (p_aOral + p_mOral)/(p_Sediment.Soil + p_Skin))

decom_tbl_long <- decom_tbl %>%
  # Pivot longer
  pivot_longer(starts_with("p_"), names_to = "Source", values_to = "Proportion")

# Turn source into a factor
decom_tbl_long$Source <- factor(decom_tbl_long$Source,
                                levels = c("p_aOral", "p_mOral", "p_Sediment.Soil", "p_Skin", "p_Unknown"))

# Order samples by oral proportion
sample_levels <- decom_tbl %>%
  arrange(Order_grouped, Species, desc(p_aOral + p_mOral)) %>% pull(new_name)

decom_tbl_long$new_name <- factor(decom_tbl_long$new_name, levels = sample_levels)
decom_tbl$new_name <- factor(decom_tbl$new_name, levels = sample_levels)

# Rearrange
decom_tbl_long <- decom_tbl_long %>% arrange(Order_grouped, new_name, Source)

# Palette for decom barplot
source_palette <- list("p_aOral"="#E54457", "p_mOral"="#AB0A1D", "p_Sediment.Soil"="#B2AD0B", "p_Skin"="#FFCC73", "p_Unknown"="#AAAAAA")

# Plot decom results
p_d <-
  ggplot(data = decom_tbl_long, aes(y = new_name, x = Proportion, fill = Source, group = Common.name)) +
  geom_bar(stat = "identity", colour = NA) +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  scale_fill_manual(values = source_palette, name = "",
                    labels = c("ancient oral", "modern oral", "sediment/soil", "skin", "unknown")) +
  scale_x_continuous(expand = c(0,0)) + 
  theme(legend.position = "top",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        plot.margin = margin(t=1, r=10, b=1, l=1)) + 
  guides(fill = guide_legend(ncol = 2))

# Plot oral to contam ratio
p_r <- ggplot(data = decom_tbl, aes(x = oral_contam_ratio, fill = oral_contam_ratio, y = new_name)) +
  geom_point(shape = 21) +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  scale_colour_viridis_c() +
  theme(legend.position="none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  # Add line at ratio = 1
  geom_vline(xintercept = 1, linetype = "dashed") +
  xlab("oral to\nsoil+skin\nratio") +
  scale_x_log10(breaks = c(0.1, 1, 10))

#### Number of reads ####

# Number of unmapped reads per sample
p_c <- ggplot(data.frame(phy_sp@sam_data), aes(x = unmapped_count, y = new_name, fill = unmapped_count)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c() +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  theme(legend.position="none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  # Add line at threshold
  geom_vline(xintercept = min_samp, linetype = "dashed") +
  xlab("# reads") +
  scale_x_log10(breaks = c(10^2, 10^5))

# Show species as a bar
species_bar <- decom_tbl %>% select(new_name, Species, Common.name, Order_grouped) %>% unique %>%
               arrange(new_name) %>%
               # Choose one label per species (in the middle of the species group)
               group_by(Common.name, Species) %>% mutate(order = row_number()) %>%
               mutate(label = ifelse(order == floor(mean(order)), Common.name, "")) %>%
               mutate(label = ifelse(is.na(label), as.character(Species), label))

p_bar <-
  ggplot(data = species_bar, aes(y = new_name, x=1, fill = Species, group = Common.name)) +
  geom_tile() + scale_fill_manual(values = species_palette, name = "") +
  facet_grid(rows = vars(Order_grouped), scales = "free", space = "free", switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank()) +
        scale_y_discrete(position = "right",
                        label = setNames(species_bar$label, species_bar$new_name)) +
    xlab("") + ylab("")

#### Run cuperdec ####
## Run cuperdec to identify good and bad samples

# Get the OTU table
otu_table <- phy_sp@otu_table %>% data.frame %>% rownames_to_column(var="Taxon")
cpdc_taxa <- load_taxa_table(otu_table) 

# Get isolation sources using the habitat table
# If it is reported as oral more frequently than environmental (soil + marine) it is considered oral
isol <- habitats %>% filter(habitat %in% c("oral", "marine", "soil")) %>%
  pivot_wider(names_from = "habitat", values_from = "occurences", id_cols = "taxon", values_fill = 0) %>%
  group_by(taxon) %>% filter(oral > (soil+marine)) %>%
  select("taxon") %>% unique %>% rename("Taxon" = "taxon") %>%
  mutate(Isolation_Source="oral")

cpdc_db <- load_database(isol, target = "oral")

# Metadata map
metadata_map <- load_map(data.frame(phy_sp@sam_data),
                         sample_col = "new_name",
                         source_col = "Order")

# Run cuperdec
curves <- calculate_curve(cpdc_taxa, database = cpdc_db)

filter_result <- simple_filter(curves, percent_threshold = 50)

cpdc_plot <- plot_cuperdec(curves, metadata_map, filter_result) + theme(plot.background=element_rect(fill="white"))

ggsave(file=file.path(subdir, "cuperdec_plot.png"), cpdc_plot, width=16, height=8)

# Add cuperdec info to main plot
cpdc_label <- data.frame(new_name = sample_levels,
                         Passed = filter_result$Passed[match(sample_levels, filter_result$Sample)]) %>%
              mutate(label = ifelse(Passed, "*", ""), oral_contam_ratio = 100) %>%
              mutate(Order_grouped = phy_sp@sam_data$Order_grouped[match(new_name, phy_sp@sam_data$new_name)])

p_r_cpdc <- p_r + geom_text(data = cpdc_label, aes(label = label), size = 3, colour = "#AB0A1D")

# Combine plots
p <- plot_grid(p_d, p_r_cpdc, p_c, p_bar, ncol = 4, align = "h", axis = "tb", rel_widths = c(1, 0.3, 0.3, 0.8))

ggsave(file=file.path(subdir, "assess_samples.png"), p, width=8, height=18)

### Combine tables
assess_samples <- decom_tbl %>% left_join(data.frame(phy_sp@sam_data) %>% select(new_name, unmapped_count), by = "new_name")
assess_samples$passed_cuperdec <- filter_result$Passed[match(assess_samples$new_name, filter_result$Sample)]

# Save table
write.table(assess_samples, file=file.path(subdir, "assess_samples.csv"), sep=",", row.names=FALSE, quote=FALSE)

###################
#### FILTERING ####
###################

# Filter out taxa than on average have a higher abundance in negatives
common_in_contols <- assess_taxa %>% filter(mean_ratio < 2) %>% pull(OTU)
phy_sp_f <- phy_sp %>% subset_taxa(!(taxa_names(phy_sp) %in% common_in_contols))

# Filter out abundances less than the threshold
phy_sp_f@otu_table[transform(phy_sp_f, "compositional")@otu_table < min_ab] <- 0

# Remove samples and taxa with no reads and negatives
phy_sp_f <- prune_samples(sample_sums(phy_sp_f) > 0 & !phy_sp_f@sam_data$is.neg, phy_sp_f)
phy_sp_f <- prune_taxa(taxa_sums(phy_sp_f) > 0, phy_sp_f)

# Get number of taxa
phy_sp_f@sam_data$taxa_filt <- estimate_richness(phy_sp_f, measures="Observed")$Observed

# Normalise
phy_sp_f_norm <- phy_sp_f %>% microbiome::transform("clr")

saveRDS(phy_sp_f, file.path(subdir, "phyloseq_objects", "phy_sp_f.RDS"))
saveRDS(phy_sp_f_norm, file.path(subdir, "phyloseq_objects", "phy_sp_f_norm.RDS"))

########################
#### CREATE SUBSETS ####
########################

## Marine terrestrial comparisons
phy_habitat <- phy_sp_f %>%
  subset_samples(Species %in% c("Otaria byronia", "Meles meles", "Ursus arctos", "Orcinus orca", "Hippopotamus amphibius", "Sus scrofa", "Dugong dugon", "Loxodonta africana"))

phy_habitat <- phy_habitat %>% subset_taxa(taxa_sums(phy_habitat) > 0)

phy_habitat_norm <- phy_habitat %>% microbiome::transform("clr")

## Planned contrasts in more deeply sampled Order (Artiodactyla, Carnivora, Primates)
# Artiodactyla
phy_artio <- phy_sp_f %>%
  subset_samples(Order == "Artiodactyla")

phy_artio <- phy_artio %>% subset_taxa(taxa_sums(phy_artio) > 0)

phy_artio_norm <- phy_artio %>% microbiome::transform("clr")

# Carnivora
phy_carni <- phy_sp_f %>%
  subset_samples(Order == "Carnivora")

phy_carni <- phy_carni %>% subset_taxa(taxa_sums(phy_carni) > 0)

phy_carni_norm <- phy_carni %>% microbiome::transform("clr")

# Primates
phy_prim <- phy_sp_f %>%
  subset_samples(Order == "Primates")

phy_prim <- phy_prim %>% subset_taxa(taxa_sums(phy_prim) > 0)

phy_norm_prim <- phy_prim %>% microbiome::transform("clr")

# All together
phy_deep <- phy_sp_f %>%
  subset_samples(Order %in% c("Artiodactyla", "Carnivora", "Primates"))

phy_deep <- phy_deep %>% subset_taxa(taxa_sums(phy_deep) > 0)

phy_norm_deep <- phy_deep %>% microbiome::transform("clr")

# Save phyloseq objects
for (obj in c("phy_habitat", "phy_artio", "phy_carni", "phy_prim", "phy_deep",
               "phy_habitat_norm", "phy_artio_norm", "phy_carni_norm", "phy_norm_prim", "phy_norm_deep")) {
  saveRDS(get(obj), file.path(subdir, "phyloseq_objects", paste0(obj, ".RDS")))
}
