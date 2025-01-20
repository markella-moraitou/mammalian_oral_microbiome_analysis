##### GENERATE PHYLOSEQ #####

#### Use MALT output table and sample metadata
#### to generate a phyloseq object for further analysis

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(renv)
library(phyloseq)
library(stringr)
library(tibble)
library(taxize)
library(microbiome)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # subdirectory for the output of this script

options(ENTREZ_KEY = Sys.getenv("API_KEY"))

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

#######################
#####  LOAD INPUT #####
#######################

#### Define file paths
# OTU table
otu_table_path <- file.path(indir, "malt_abundance_matrix_sam.txt")
tax_names = file.path(indir, "taxname_to_id.txt")

# Sample metadata
metadata_path <- file.path(indir, "sample_metadata.csv")
sample_tax_path <- file.path(indir, "host_taxonomy.csv")
species_habitats_path <- file.path(indir, "species_habitats.csv")
quant_diet_path <- file.path(indir, "Lintulaakso_diet_filtered.csv")
rc_path <- file.path(indir, "read_count.csv")
rl_path <- file.path(indir, "read_length.csv")
decom_path <- file.path(indir, "decOM_output.csv")

# HOMD database
#homd_path <- file.path(projdir, "input", "HOMD_taxon_table.csv")

#### Load files
# Load OTU table
otu_table <- read.table(otu_table_path, sep="\t", comment.char="", header=TRUE)
names_to_ids <- read.table(tax_names, sep="\t", comment.char="", header=FALSE, col.names=c("names", "ids"), colClasses=c("character", "character"))

# Load metadata
metadata <- read.csv(metadata_path) # Sample metadata
sample_tax <- read.csv(sample_tax_path) %>%  # Sample taxonomy
  select(museum.species, genus, subfamily, infraorder, suborder, order, superorder, Common.name)
species_habitats <- read.csv(species_habitats_path) # Species habitats
quant_diet <- read.csv(quant_diet_path) # Diet quantification data from Lintulaakso et al. 2023 paper
rc <- read.csv(rc_path) %>% rename_with( ~ paste0(., "_count")) # read count per step
rl <- read.csv(rl_path) %>% rename_with( ~ paste0(., "_avlength")) # average read length per step
decom <- read.csv(decom_path, row.names = NULL) %>% # decOM output
  # remove entries where no kmers have been counted
  filter(rowSums(!is.na(select(., starts_with("p_")))) > 0)

# Load HOMD database
#homd <- read.csv(homd_path, skip=1) %>% filter(Body_site == "Oral") # Keep only oral taxa

###############################################
#### COMBINE AND TIDY UP SAMPLE METADATA  #####
###############################################

# Normalise colnames
colnames(sample_tax) <- str_to_title(colnames(sample_tax))

# Combine all sample and host species metadata in one big table
meta <- metadata

meta <-
  meta %>%
  left_join(sample_tax, relationship = "many-to-many", by=c("Species"="Museum.species")) %>%
  left_join(species_habitats, relationship = "many-to-many") %>%
  left_join(quant_diet, relationship="many-to-many") %>%
  left_join(rc, by=c("Ext.ID"="sample_count")) %>%
  left_join(rl, by=c("Ext.ID"="sample_avlength")) %>% 
  left_join(decom, by=c("Ext.ID"="Sink")) %>%
  # Combine different rows for pooled samples
  group_by(Ext.ID) %>%
  summarise_all(.funs = function(.cols) {
    if (n_distinct(.cols) > 1) {
      paste(.cols, collapse = "&")
    } else {
      .cols[1]
    }
  }) %>%
  rename("diet.general"="calculated_cluster_main_diet") %>%
  as.data.frame

# Turn order and species to a factor and create a separate order for blanks and controls
meta$Order <- ifelse(meta$Species %in% c("Environmental control", "Extraction blank", "Library blank"),
                    "Control/blank", meta$Order)

ord_levels <- unique(sort(meta$Order))
ord_levels <- ord_levels[ord_levels != c("Control/blank")] %>% append("Control/blank")
meta$Order <- factor(meta$Order, levels=ord_levels)

spe_levels <- meta %>% arrange(Order, Species) %>% pull(Species) %>% unique
meta$Species <- factor(meta$Species, levels=spe_levels)

# Create a new Order column with rare orders grouped together for easier plotting
meta <- meta %>% group_by(Order) %>%
  mutate(Order_grouped = case_when(n_distinct(Species) > 1 ~ Order,
                                   TRUE ~ "Rest")) %>% as.data.frame

levels <- unique(meta$Order_grouped)
levels <- levels[which(! levels %in% c("Rest", "Control/blank"))] %>% append(c("Rest", "Control/blank"))

meta$Order_grouped <- factor(meta$Order_grouped, levels = levels)

## Get better samples names
rename <- meta %>%
  # Give common name to blanks and controls
  mutate(Common.name = case_when(Species %in% c("Environmental control", "Extraction blank", "Library blank") ~ Species, TRUE ~ Common.name)) %>%
  select(Ext.ID, Common.name) %>% rename("old_name"="Ext.ID") %>%
  # new names will consist of the first letter of the genus, the first three of the species epithet and a number
  group_by(Common.name) %>% mutate(num=row_number() %>% str_pad(width = 2, pad = "0")) %>%
  separate(col=Common.name, into=c("part1", "part2", "part3", "part4"), fill="left", sep=" ") %>%
  # Turn NAs to ""
  mutate(part1=ifelse(is.na(part1), "", part1),
         part2=ifelse(is.na(part2), "", part2),
         part3=ifelse(is.na(part3), "", part3),
         part4=ifelse(is.na(part4), "", part4)) %>%
  # Use first 4 letters of last adjective (part3) and the entire last word (part4)
  mutate(new_name=paste(str_to_lower(ifelse(str_length(part3)> 3, str_sub(part3, 1, 1), part3)),
                        str_to_title(part4), "_", num, sep="")) %>%
  select(old_name, new_name)

# Temporary bit for this unknown library
rename$new_name[rename$old_name=="maybe_DM_017"] <- paste0(rename$new_name[rename$old_name == "DM_017"], ".maybe")

meta$new_name <- rename$new_name[match(meta$Ext.ID, rename$old_name)]

# Keep a version with all lab metadata (not just for samples in OTU table)
meta_all <- meta

############################
#### PREPARE OTU TABLE #####
############################

#### Adapt table to use in phyloseq ####
# Get taxids (this will be used to get taxonomic info)
tbl <- otu_table
rownames(tbl) <- names_to_ids$ids[match(rownames(tbl), names_to_ids$names)]

# Get new names
tbl <- tbl %>% rename_with(~str_remove(., "_unmapped"), everything()) %>%
  rename_with(~rename$new_name[match(., rename$old_name)], everything())

# Keep only samples that exist in tbl
meta <- meta %>% filter(new_name %in% colnames(tbl))

# Get rownames
rownames(meta) <- meta$new_name

###############################
####  GET TAXONOMIC INFO  #####
###############################

# If the file already exists, load it, otherwise start an empty file
if(file.exists(file.path(subdir, "taxonomy_all.csv"))) {
  taxonomy <- read.csv(file.path(subdir, "taxonomy_all.csv"), colClasses = "character")
  cont <- TRUE
} else {
  taxonomy <- data.frame(superkingdom=character(),
                       phylum=character(),
                       order=character(),
                       family=character(),
                       genus=character(),
                       species=character(),
                       tax.id=character())
  cont <- FALSE
}

# Retrieve taxonomic information from NCBI
taxids <- rownames(tbl) # Tax IDs from OTU table
levels <- colnames(taxonomy) # Taxonomic levels

for (i in 1:length(taxids)) {
  id <- taxids[i]
  # If a pre-existing file was loaded check if the id is already there and if yes skip
  if (cont==TRUE) {
    if (id %in% taxonomy$tax.id) {
      cat(i, "ID ", id, " already there, skipping\n")
      next
    }
  }
  # Get classification of this taxon.
  #If it fails, retry after 5 seconds, if it fails a second time, return NA
  classif <- tryCatch({
    taxize::classification(id, db="ncbi")
  }, error = function(e) {
    cat("Error occurred for taxon ", id, ": \n", conditionMessage(e), "\n")
    cat("Retrying after 5 seconds...\n")
    Sys.sleep(5)
    tryCatch({
      taxize::classification(id, db="ncbi")
    }, error = function(e) {
      cat("Error occurred again for taxon ", id, ": \n", conditionMessage(e), "\n")
      cat("Returning NA \n")
      na_table = data.frame(names = NA, rank = levels, name = NA)
      return(na_table)
    })
  })
  # Bring to the right format to combine with the big taxonomy table
  classif_tbl <- data.frame()
  for (level in levels){
    # Get the name on that taxonomic level
    name <- classif[[1]] %>% filter(rank == level) %>% pull(name)
    # If there is no name at that level, just use NA
    name <- ifelse(is.null(name), as.character(NA), name)
    classif_tbl[1, level] <- name
  }
  classif_tbl["tax.id"] <- id
  
  # Add to the big taxonomy table
  taxonomy <- rbind(taxonomy, classif_tbl)
  
  # Print most specific taxonomic info as it progresses
  cat(i, ":", length(taxids), "ID ", id, " -- ", classif_tbl$species, "\n")
}

taxonomy <- as.matrix(taxonomy)
# Get taxa names as row names
row.names(taxonomy) <- taxonomy[, "tax.id"]

###################################
####  CREATE PHYLOSEQ TABLES  #####
###################################

#### Create phyloseq object, agglomerate at species level, and filter
OTU = otu_table(tbl, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy) %>% unique
SAM = sample_data(meta)

# Create phyloseq object and filter out empty samples
cat("Creating phyloseq table\n")
phy <- phyloseq(OTU, SAM, TAX)

phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- subset_samples(phy, sample_sums(phy) > 0)

# Agglomerate to species level and keep only prokaryotes and archaea
cat("Agglomerating to species level\n")
phy_sp <- phy %>% tax_glom(taxrank="species") %>%
 subset_taxa(superkingdom %in% c("Bacteria", "Archaea"))

phy_sp <- subset_taxa(phy_sp, taxa_sums(phy_sp) > 0)
phy_sp <- subset_samples(phy_sp, sample_sums(phy_sp) > 0)

# Get taxa names instead of ids
taxa_names(phy_sp) <- as.vector(phy_sp@tax_table[,"species"])

#### Collect some metadata

# Collect number of OTUs per sample
phy_sp@sam_data$taxa_raw <- estimate_richness(phy_sp, measures="Observed")$Observed

#Add column indicating samples and controls
phy_sp@sam_data$is.neg <- grepl("blank|control", phy_sp@sam_data$Order_grouped)

# Calculate oral to soil ratio according to DecOM results
phy_sp@sam_data <- phy_sp@sam_data %>% data.frame %>% mutate(oral_to_soil_ratio=(p_mOral + p_aOral)/p_Sediment.Soil) %>% sample_data

# CLR-normalisation
phy_sp_norm <- phy_sp %>% transform('clr')

#####################
#### SAVE OUTPUT ####
#####################

# Phyloseq objects
dir.create(file.path(subdir, "phyloseq_objects"), showWarnings = FALSE)

saveRDS(phy_sp, file.path(subdir, "phyloseq_objects", "phy_sp.RDS"))
saveRDS(phy_sp_norm, file.path(subdir, "phyloseq_objects", "phy_sp_norm.RDS"))

# Constituent parts of phyloseq object
write.table(otu_table(phy_sp), file.path(subdir, "phyloseq_objects", "phy_sp_OTU.csv"), sep=",", row.names=TRUE, quote=FALSE)
write.table(tax_table(phy_sp), file.path(subdir, "phyloseq_objects", "phy_sp_TAX.csv"), sep=",", row.names=TRUE, quote=FALSE)
write.table(data.frame(sample_data(phy_sp)), file.path(subdir, "phyloseq_objects", "phy_sp_SAM.csv"), sep=",", row.names=TRUE, quote=FALSE)
write.table(otu_table(phy_sp_norm), file.path(subdir, "phyloseq_objects", "phy_sp_norm_OTU.csv"), sep=",", row.names=TRUE, quote=FALSE)

# Entire taxonomy table
write.table(taxonomy, file.path(subdir, "taxonomy_all.csv"), sep=",", row.names=FALSE, quote=TRUE)

# Entire metadata table
# Get a big metadata table with all lab samples
meta_all <- meta_all %>% left_join(data.frame(phy_sp@sam_data)) %>%
  mutate(is.neg=grepl("blank|control", Species))

write.table(meta_all, file.path(subdir, "metadata_all.csv"), sep=",", row.names=FALSE, quote=FALSE)

# Save names and ids of retained taxa
names_to_ids_filt <- names_to_ids %>% filter(names %in% rownames(otu_table(phy_sp)))
write.table(names_to_ids_filt, file.path(subdir, "names_to_ids_filt.csv"), sep=",", row.names=FALSE, quote=FALSE)
