#### FUNCTIONS TO USE WITH PCA ORDINATIONS ####


## Add species centroid as phylopics
centroids <- function(ordination, phyloseq) {
  # Get ordination vectors as data frame and species info
  ord_df <- data.frame(scores(ordination, choices = c(1:4), display = "sites"))
  ord_df <- ord_df %>%
    cbind(select(data.frame(phyloseq@sam_data), c(Common.name, Species, Order_grouped, Order, diet.general, habitat.general)))
  # Group and calculate means
  centroids <- ord_df %>% group_by(Common.name, Species, Order_grouped, Order, diet.general, habitat.general) %>%
          summarise_all(.funs = mean, na.rm = TRUE)
  centroids <- left_join(centroids, phylopics, by = c("Species" = "Species"))
  return(centroids)
}

## Get loadings plot to display alongside PCA
loadings_plot <- function(ordination, axes, top_taxa = 20) {
  # Get loadings
  loadings <- data.frame(scores(ordination, choices = axes, display = "species"))
  axes_names=paste0("PC", axes)
  eigval <- ordination$CA$eig[axes]
  # Keep only taxa with highest loadings
  loadings_filt <- loadings %>%
    # Arrange by the sum of the absolute values of the loadings, weighted by the eigenvalues (longest arrow)
    arrange(desc(abs(!!sym(axes_names[1])*eigval[1] + !!sym(axes_names[2])*eigval[2]))) %>%
    # Keep only top taxa
    head(top_taxa) %>% rownames_to_column("taxon") %>%
    # Get labels for plotting
    mutate(genus = gsub(" .*", "", taxon)) %>%
    # Create labels for taxa by keeping only first letter of genus
    mutate(label = gsub("([A-Z])[a-z]+", "\\1.", taxon)) %>%
    # Arrange by first axis
    arrange(desc(!!sym(axes_names[1]))) %>% mutate(label = factor(label, levels = label))
  # Keep top 8 genera and group the rest as "Other"
  top_genera <- table(loadings_filt$genus) %>% sort(decreasing = TRUE) %>% names() %>% head(8)
  loadings_filt <- loadings_filt %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other"))
  # Get genus palette
  genera <- loadings_filt %>% filter(genus != "Other") %>% pull(genus) %>% unique
  genus_palette <- setNames(wes_palette("Darjeeling1", length(genera), type = "continuous"), genera)
  genus_palette["Other"] <- "grey"
  # Set genus as factor
  loadings_filt$genus <- factor(loadings_filt$genus, levels = names(genus_palette))
  # Melt and plot
  loadings_melt <- pivot_longer(loadings_filt, cols = starts_with("PC"), names_to = "PC", values_to = "loading")
  p <- ggplot(loadings_melt, aes(y = label, x = loading, fill = genus)) +
        geom_bar(stat = "identity", colour = "black") +
        facet_grid(~PC, scales = "free_y", space = "free_y", switch = "y") +
        geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
        scale_fill_manual(values = genus_palette, name = "Genus") +
        theme(axis.text.y = element_text(face = "italic", family = "serif", size = 8),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.grid = element_blank(),
              strip.background = element_rect(colour = "white"),
              legend.position = "left")
  return(p)
}

## Get loadings plot to display alongside RDA
loadings_plot_rda <- function(ordination, axes, top_taxa = 20) {
  # Get loadings
  loadings <- data.frame(scores(ordination, choices = axes, display = "species"))
  axes_names=paste0("RDA", axes)
  eigval <- ordination$CA$eig[axes]
  # Keep only taxa with highest loadings
  loadings_filt <- loadings %>%
    # Arrange by the sum of the absolute values of the loadings, weighted by the eigenvalues (longest arrow)
    arrange(desc(abs(!!sym(axes_names[1])*eigval[1] + !!sym(axes_names[2])*eigval[2]))) %>%
    # Keep only top taxa
    head(top_taxa) %>% rownames_to_column("taxon") %>%
    # Get labels for plotting
    mutate(genus = gsub(" .*", "", taxon)) %>%
    # Create labels for taxa by keeping only first letter of genus
    mutate(label = gsub("([A-Z])[a-z]+", "\\1.", taxon)) %>%
    # Arrange by first axis
    arrange(desc(!!sym(axes_names[1]))) %>% mutate(label = factor(label, levels = label))
  # Keep top 8 genera and group the rest as "Other"
  top_genera <- table(loadings_filt$genus) %>% sort(decreasing = TRUE) %>% names() %>% head(8)
  loadings_filt <- loadings_filt %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other"))
  # Get genus palette
  genera <- loadings_filt %>% filter(genus != "Other") %>% pull(genus) %>% unique
  genus_palette <- setNames(wes_palette("Darjeeling1", length(genera), type = "continuous"), genera)
  genus_palette["Other"] <- "grey"
  # Set genus as factor
  loadings_filt$genus <- factor(loadings_filt$genus, levels = names(genus_palette))
  # Melt and plot
  loadings_melt <- pivot_longer(loadings_filt, cols = starts_with("RDA"), names_to = "RDA", values_to = "loading")
  p <- ggplot(loadings_melt, aes(y = label, x = loading, fill = genus)) +
        geom_bar(stat = "identity", colour = "black") +
        facet_grid(~RDA, scales = "free_y", space = "free_y", switch = "y") +
        geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
        scale_fill_manual(values = genus_palette, name = "Genus") +
        theme(axis.text.y = element_text(face = "italic", family = "serif", size = 8),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.grid = element_blank(),
              strip.background = element_rect(colour = "white"),
              legend.position = "left")
  return(p)
}
