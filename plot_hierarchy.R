plot_hierarchy <- function(physeq, number_of_taxa_to_plot=50, plot_type="sankey"){
  
  library(phyloseq)
  library(tidyverse)
  library(hiervis)
  
  # start by defining abundances and top taxa functions. borrowed from microbiome package https://doi.org/doi:10.18129/B9.bioc.microbiome
  abundances <- function (x, transform = "identity") 
  {
    if (class(x) == "phyloseq") {
      otu <- as(otu_table(x), "matrix")
      if (!taxa_are_rows(x) && ntaxa(x) > 1 && nsamples(x) > 
          1) {
        otu <- t(otu)
      }
      if (ntaxa(x) == 1) {
        otu <- matrix(otu, nrow = 1)
        rownames(otu) <- taxa(x)
        colnames(otu) <- sample_names(x)
      }
      if (nsamples(x) == 1) {
        otu <- matrix(otu, ncol = 1)
        rownames(otu) <- taxa(x)
        colnames(otu) <- sample_names(x)
      }
    }
    else if (is.vector(x)) {
      otu <- as.matrix(x, ncol = 1)
    }
    else {
      otu <- x
    }
    if (!transform == "identity") {
      otu <- transform(otu, transform)
    }
    otu
  }
  
  
  top_taxa <- function (x, n = ntaxa(x)) 
  {
    names(sort(rowSums(abundances(x)), decreasing = TRUE)[seq_len(n)])
  }
  
  # first, we prune to the specified number of taxa
  physeq_abund <-  prune_taxa(top_taxa(physeq, number_of_taxa_to_plot), physeq) 
  
  # then convert Kingdom-Genus tax table to dataframe
  tax_df <- as.data.frame(phyloseq::tax_table(physeq_abund)[,1:6])
  
  tax_df$Life_type <- rep("Cellular", nrow(tax_df))
  
  # add "path" column for each genus; this is used by hiervis to generate plot
  tax_df$path <- paste(tax_df$Life_type, tax_df$Kingdom, tax_df$Phylum, tax_df$Class, tax_df$Order, tax_df$Family, tax_df$Genus, sep = "/") 
  
  #populate data frame with taxa abundances
  tax_df$Abundance <- as.numeric(taxa_sums(physeq_abund))
  
  # subset to only include path and Abundance columns
  tax_df_subsetted <- tax_df[,c("path", "Abundance")]
  
  # "path" column also requires this hierarchical information
  life_type <- "Cellular"
  k <- unique(paste(tax_df$Life_type, tax_df$Kingdom, sep = "/"))
  kp <- unique(paste(tax_df$Life_type,tax_df$Kingdom, tax_df$Phylum,  sep = "/"))
  kpc <- unique(paste(tax_df$Life_type,tax_df$Kingdom, tax_df$Phylum, tax_df$Class,  sep = "/"))
  kpco  <- unique(paste(tax_df$Life_type,tax_df$Kingdom, tax_df$Phylum, tax_df$Class, tax_df$Order,  sep = "/"))
  kpcof  <- unique(paste(tax_df$Life_type,tax_df$Kingdom, tax_df$Phylum, tax_df$Class, tax_df$Order, tax_df$Family, sep = "/"))
  kpcofg   <- unique(paste(tax_df$Life_type,tax_df$Kingdom, tax_df$Phylum, tax_df$Class, tax_df$Order, tax_df$Family, tax_df$Genus, sep = "/"))
  
  #combines hierarchical info
  path_header <- c(life_type, k, kp, kpc, kpco, kpcof, kpcofg)
  
  # makes data frame from hierarchical info
  path_header_df <- data.frame(path_header, rep(NA, length(path_header)))
  colnames(path_header_df) <- c("path", "Abundance")
  
  # combine hierarchical info and abundance data frames
  df_to_plot <- rbind(path_header_df, tax_df_subsetted)
  
  return(hiervis( df_to_plot,  plot_type, nameField = "path", pathSep = "/", valueField = "Abundance", stat = "sum"))
}
