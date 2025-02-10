##### NORMALISATIONS #####

#### Access OTU tree from Open Tree of Life, perform philr and clr normalisations
#### and create subsets for further analysis

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(tidyr)
library(rotl)
library(ape)
library(stringr)
library(philr)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it does not exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))

#######################
#### GET TAXA TREE ####
#######################

## Match phyloseq taxa to Open Tree of Life 

# For a few taxa, we need to use a synonym to match
synonyms <- list("Hoylesella enoeca" = "Prevotella enoeca",
                 "Segatella bryantii" = "Prevotella bryantii",
                 "Segatella oris" = "Prevotella oris",
                 "Xylanibacter ruminicola" = "Prevotella ruminicola",
                 "Flagellimonas marinaquae" = "Muricauda aquimarina",
                 "Thermaerobacillus caldiproteolyticus" = "Anoxybacillus caldiproteolyticus",
                 "Aedoeadaptatus ivorii" = "Peptoniphilus ivorii")

synonyms_df <- data.frame(phy_name = names(synonyms), search_string = unlist(synonyms))

write.table(synonyms_df, file=file.path(outdir, "rotl_synonyms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Match names to open tree taxonomic names
bacteria_names <- phy_sp_f %>% subset_taxa(superkingdom == "Bacteria") %>% taxa_names
# Replace with synonyms
bacteria_names <- ifelse(bacteria_names %in% synonyms_df$phy_name, synonyms_df$search_string[match(bacteria_names, synonyms_df$phy_name)], bacteria_names)

resolved_names_b <- tnrs_match_names(bacteria_names, context_name = "Bacteria")
resolved_names_b$superkingdom <- "Bacteria"

archaea_names <- phy_sp_f %>% subset_taxa(superkingdom == "Archaea") %>% taxa_names
resolved_names_a <- tnrs_match_names(archaea_names, context_name = "Archaea")
resolved_names_a$superkingdom <- "Archaea"

# Combine the bacteria and archea tables
resolved_names <- rbind(resolved_names_b, resolved_names_a)

# Remove incertae sedis and replace with genus name (better to get the genus level placement than nothing)
incertae_sedis <- resolved_names %>% filter(grepl("incertae_sedis", flags) | grepl("unplaced", flags))
incertae_sedis$genus <- gsub(" .*", "", incertae_sedis$search_string)

write.table(incertae_sedis, file = file.path(outdir, "rotl_incertae_sedis.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)

for (i in 1:nrow(incertae_sedis)) {
  search_string <- incertae_sedis$search_string[i]
  genus <- incertae_sedis$genus[i]
  superkingdom <- incertae_sedis$superkingdom[i]
  df <- tnrs_match_names(genus, context_name = superkingdom)
  df$superkingdom <- superkingdom
  resolved_names <- rbind(resolved_names, df)
}

resolved_names <- resolved_names %>% filter(!(search_string %in% incertae_sedis$search_string)) %>%
   filter(!grepl("unplaced", flags))

write.table(resolved_names, file.path(phydir, "rotl_resolved_names.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Get the phylogenetic tree
tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)

#########################
#### MANIPULATE TREE ####
#########################

## So that the tips match the phyloseq object

rename_tip_labels <- data.frame(original =  tree$tip.label,
                                ott = as.numeric(gsub(".*ott", "", tree$tip.label))) %>%
                      # Get the taxa names by matching the OTT
                      left_join(select(resolved_names, c(search_string, ott_id)), by = c("ott" = "ott_id")) %>%
                      mutate(taxon = str_to_sentence(search_string)) %>%
                      # Get only genus and species epithet (remove infraspecific names)
                      separate(taxon, c("genus", "species"), sep = " ", extra = "drop") %>%
                      mutate(species = ifelse(is.na(species), "", species),
                             genus = ifelse(is.na(genus), "", genus)) %>%
                      mutate(taxon = paste(genus, species, sep = " ")) %>%
                      # Change back to match the phyloseq taxa names
                      mutate(taxon = case_when(taxon %in% synonyms_df$search_string ~ synonyms_df$phy_name[match(taxon, synonyms_df$search_string)],
                                               TRUE ~ taxon)) %>%
                      # Check if the taxon names match those in the phyloseq object
                      mutate(match = taxon %in% taxa_names(phy_sp_f))

# Which taxa don't match those in phyloseq
no_match_tips <- rename_tip_labels %>% filter(!match) %>% select(original, ott, search_string, taxon)

# Which phyloseq taxa are not represented?
unrepr_taxa <- setdiff(taxa_names(phy_sp_f), rename_tip_labels$taxon)

for (i in 1:nrow(no_match_tips)) {
    # If there is no taxon name, try to get the full name using the OTT
    if (no_match_tips$taxon[i] == " ") {
      full_name <- tol_node_info(no_match_tips$ott[i])$taxon$unique_name
      name_split <- str_split(full_name, pattern = " ")[[1]]
      taxon_name <- paste(name_split[1], name_split[2], sep = " ") %>% str_remove(" NA")
    } else {
      taxon_name <- no_match_tips$taxon[i]
    }
    # Check if the taxon name is in the list of unrepresented taxa
    if (!(str_remove(taxon_name, " sp.") %in% unrepr_taxa)) {
      cat("Couldn't find a match for ", taxon_name, ". Checking genus...\n", sep="")
      genus_match <- unrepr_taxa[grepl(str_remove(taxon_name, pattern = " sp."), unrepr_taxa)]
      if (length(genus_match) == 1) {
        taxon_name <- genus_match
        cat("Match found: ", genus_match, "\n", sep = "")
      } else {
        cat("Error: ", length(genus_match), " matches found\n", sep = "")
        break      
      }
    }
    no_match_tips$taxon[i] <- taxon_name
}

write.table(no_match_tips, file = file.path(outdir, "rotl_no_match_tips.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Add the new names to the rename_tip_labels table
rename_tip_labels <- rename_tip_labels %>%
                      mutate(taxon = case_when(!match ~ no_match_tips$taxon[match(ott, no_match_tips$ott)],
                                              TRUE ~ taxon)) %>%
                      # Check if the taxon names match those in the phyloseq object
                      mutate(match = taxon %in% taxa_names(phy_sp_f))

write.table(select(rename_tip_labels, c(original, ott, taxon)), file = file.path(outdir, "rotl_rename_tips.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## Also rename the node labels, for better readability
# I will replace the otu codes with the actual taxon names
rename_nodes <- data.frame(original =  tree$node.label) %>%
                mutate(temp = str_remove(original, "mrcaott")) %>%
                mutate(temp = gsub(" ott", "ott", temp)) %>%
                mutate(temp = gsub(" \\(.*?\\)", "", temp)) %>%
                separate(temp, c("ott1", "ott2"), sep = "ott")

rename_nodes$new_name <- NA
for (i in 1:nrow(rename_nodes)) {
    cat("Processing node ", i, " of ", nrow(rename_nodes), "\n", sep = "")
    # For the ott1 and ott2 columns separately, check if an OTT is there, instead of a taxon name,
    # and if so, retrieve the taxon name
    if (!grepl("[a-zA-Z]", rename_nodes$ott1[i])) {
        taxon1 <- tol_node_info(rename_nodes$ott1[i])$taxon$unique_name
    } else { taxon1 <- rename_nodes$ott1[i] }
    if (!grepl("[a-zA-Z]", rename_nodes$ott2[i])) {
        taxon2 <- tol_node_info(rename_nodes$ott2[i])$taxon$unique_name
    } else { taxon2 <- rename_nodes$ott2[i] }
    rename_nodes$new_name[i] <- paste(taxon1, taxon2, sep = "-")
}

rename_nodes$new.name <- make.unique(rename_nodes$new_name)
write.table(rename_nodes, file = file.path(outdir, "rotl_rename_nodes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Relabel the tree and add to phyloseq
tree$tip.label <- rename_tip_labels$taxon[match(tree$tip.label, rename_tip_labels$original)]
tree$node.label <- rename_nodes$new_name[match(tree$node.label, rename_nodes$original)]

write.tree(tree, file = file.path(phydir, "phy_tree.tree"))

phy_sp_f <- phyloseq(otu_table(phy_sp_f), sample_data(phy_sp_f), tax_table(phy_sp_f), phy_tree(tree))

# CLR transformation
phy_sp_f_clr <- phy_sp_f %>% microbiome::transform("clr")

# PhiLR transformation
philr_otu <- philr(phy_sp_f, pseudocount=10^-5)

write.table(philr_otu, file = file.path(phydir, "phy_sp_philr_OTU.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)

phy_sp_philr <- phyloseq(otu_table(philr_otu, taxa_are_rows = FALSE), sample_data(phy_sp_f))

# Save phyloseq objects
saveRDS(phy_sp_f, file.path(phydir, "phy_sp_f.RDS"))
saveRDS(phy_sp_f_clr, file.path(phydir, "phy_sp_f_clr.RDS"))
saveRDS(phy_sp_philr, file.path(phydir, "phy_sp_philr.RDS"))

########################
#### CREATE SUBSETS ####
########################

# Create function to mirror subset of phyloseq object on other transformed objects
mirror_subset <- function(phy_subset, phy_norm, taxa = FALSE) {
  phy_norm_new <- phy_norm %>% prune_samples(sample_names(phy_norm) %in% sample_names(phy_subset), .)
  if (taxa) {
    cat("Also subsetting taxa\n")
    phy_norm_new <- phy_norm_new %>% prune_taxa(taxa_names(phy_norm) %in% taxa_names(phy_subset), .)
  }
  return(phy_norm_new)
}

## Wild animals only
phy_wild <- phy_sp_f %>% subset_samples(!(Species %in% c("Ovis aries", "Equus caballus", "Sus scrofa domesticus")))
phy_wild <- phy_wild %>% subset_taxa(taxa_sums(phy_wild) > 0)

phy_wild_clr <- mirror_subset(phy_subset=phy_wild, phy_norm=phy_sp_f_clr, TRUE)
phy_wild_philr <- mirror_subset(phy_subset=phy_wild, phy_norm=phy_sp_philr, FALSE)

## Marine terrestrial comparisons
phy_habitat <- phy_sp_f %>%
  subset_samples(Species %in% c("Otaria byronia", "Meles meles", "Ursus arctos", "Orcinus orca", "Hippopotamus amphibius", "Sus scrofa", "Dugong dugon", "Loxodonta africana"))
phy_habitat <- phy_habitat %>% subset_taxa(taxa_sums(phy_habitat) > 0)

phy_habitat_clr <- mirror_subset(phy_subset=phy_habitat, phy_norm=phy_sp_f_clr, TRUE)
phy_habitat_philr <- mirror_subset(phy_subset=phy_habitat, phy_norm=phy_sp_philr, FALSE)

## Planned contrasts in more deeply sampled Order (Artiodactyla, Carnivora, Primates)
# Artiodactyla
phy_artio <- phy_sp_f %>%
  subset_samples(Order == "Artiodactyla")
phy_artio <- phy_artio %>% subset_taxa(taxa_sums(phy_artio) > 0)

phy_artio_clr <- mirror_subset(phy_subset=phy_artio, phy_norm=phy_sp_f_clr, TRUE)
phy_artio_philr <- mirror_subset(phy_subset=phy_artio, phy_norm=phy_sp_philr, FALSE)

# Carnivora
phy_carni <- phy_sp_f %>%
  subset_samples(Order == "Carnivora")
phy_carni <- phy_carni %>% subset_taxa(taxa_sums(phy_carni) > 0)

phy_carni_clr <- mirror_subset(phy_subset=phy_carni, phy_norm=phy_sp_f_clr, TRUE)
phy_carni_philr <- mirror_subset(phy_subset=phy_carni, phy_norm=phy_sp_philr, FALSE)

# Primates
phy_prim <- phy_sp_f %>%
  subset_samples(Order == "Primates")
phy_prim <- phy_prim %>% subset_taxa(taxa_sums(phy_prim) > 0)

phy_prim_clr <- mirror_subset(phy_subset=phy_prim, phy_norm=phy_sp_f_clr, TRUE)
phy_prim_philr <- mirror_subset(phy_subset=phy_prim, phy_norm=phy_sp_philr, FALSE)

# All together
phy_deep <- phy_sp_f %>%
  subset_samples(Order %in% c("Artiodactyla", "Carnivora", "Primates"))
phy_deep <- phy_deep %>% subset_taxa(taxa_sums(phy_deep) > 0)

phy_deep_clr <- mirror_subset(phy_subset=phy_deep, phy_norm=phy_sp_f_clr, TRUE)
phy_deep_philr <- mirror_subset(phy_subset=phy_deep, phy_norm=phy_sp_philr, FALSE)

# Save phyloseq objects
for (obj in c("phy_wild", "phy_habitat", "phy_artio", "phy_carni", "phy_prim", "phy_deep",
               "phy_wild_clr", "phy_habitat_clr", "phy_artio_clr", "phy_carni_clr", "phy_prim_clr", "phy_deep_clr",
               "phy_wild_philr", "phy_habitat_philr", "phy_artio_philr", "phy_carni_philr", "phy_prim_philr", "phy_deep_philr")) {
  saveRDS(get(obj), file.path(phydir, paste0(obj, ".RDS")))
}
