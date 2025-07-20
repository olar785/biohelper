#' Loads normalized_otu_taxo
#'
#' @description
#' This function normalises taxonomic assignments of ASVs belonging to the same OTU by using
#' the taxonomic assignment of the ASV that has the highest taxonomic resolution.
#'
#'
#' @param
#' physeq_obj = A phyloseq object
#' @param
#' otu_column = Name of the OTU column
#' @param
#' asv_column = Name of the ASV column
#'
#' @export
#' @examples
#' # Normalizing taxonomy followed by merging ASVs to OTUs
#' normalized_otu_taxo(physeq_obj = ps_obj) %>% speedyseq::tax_glom(taxrank = "otu")
#'
#'
normalized_otu_taxo <- function(physeq_obj, otu_column = "otu", asv_column = "asv") {
  # Extract taxonomic table as a data frame
  tax_table_df <- as.data.frame(phyloseq::tax_table(physeq_obj))

  # Function to find the ASV with the highest taxonomic resolution within each OTU
  get_highest_resolution_taxa <- function(subset_df) {
    # Count the number of non-NA entries for each ASV
    num_non_na <- rowSums(!is.na(subset_df))

    # Get the ASV with the maximum number of non-NA taxonomic levels
    highest_resolution_asv <- subset_df[which.max(num_non_na), ]

    # Return the taxonomy of the highest resolution ASV
    return(highest_resolution_asv)
  }

  # Group taxa by the OTU level
  tax_table_by_otu <- tax_table_df %>%
    dplyr::group_by_at(otu_column) %>%
    dplyr::summarise_all(~ if (all(is.na(.))) NA else .[which.max(!is.na(.))])

  # Get unique OTUs
  unique_otus <- unique(tax_table_by_otu[[otu_column]])

  # Initialize a progress bar
  pb <- txtProgressBar(min = 0, max = length(unique_otus), style = 3)

  # Loop over each OTU and update the progress bar
  for (i in seq_along(unique_otus)) {
    otu <- unique_otus[i]

    # Get the taxa with the highest resolution for this OTU
    highest_taxa <- get_highest_resolution_taxa(tax_table_df[tax_table_df[[otu_column]] == otu, ])

    # Assign the highest resolution taxonomy to all ASVs that belong to this OTU
    tax_table_df[tax_table_df[[otu_column]] == otu, !(colnames(tax_table_df) %in% "asv")] <- highest_taxa %>% dplyr::select(-asv)

    # Update the progress bar
    setTxtProgressBar(pb, i)
  }

  # Close the progress bar
  close(pb)

  # Return the modified phyloseq object with updated taxonomy table
  physeq_obj@tax_table <- phyloseq::tax_table(as.matrix(tax_table_df))

  return(physeq_obj)
}

