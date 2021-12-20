#' Loads simple_decon
#'
#' This function takes in a phyloseq object with sample data and performs a simple decontamination of ASV data
#' by either completely removing contaminant ASVs and their associated reads across the dataset ("complete_asv_removal")
#' or by subtracting read count per ASV using their max read count in blank(s) ("max_v"). The latter option does not
#' take into account the compositional nature of the data. If using the latter option,
#' the user may want to perform this step per batch of sample (e.g. DNA extraction or PCR batch).
#' The sample_data MUST have a column labeled sample_id and a column labeled sample_type.
#' Non blank samples must be labeled "sample".
#'For example:
#'sample_id   blank_type            DNA_extraction_batch      extraction_method   etc.
#'sample1     sample                1                         manual
#'sample2     sampling_blank        1                         manual
#'sample3     dna_extraction_blank  2                         robot
#'sample4     sample                2                         robot
#'sample5     pcr_blank             2                         robot
#'
#' @export
#' @examples
#' simple_decon(ps_test_data, method = "complete_asv_removal")
#' simple_decon(ps_test_data, method = "max_v")

simple_decon= function(ps_obj, method = "complete_asv_removal"){
  #ps_obj@sam_data$sample_type = suppressWarnings(str_replace_all(ps_obj@sam_data$sample_type, pattern = c("^NA$","^na$","^Na$","^NaN$","^nan$",""), replacement = NA_character_))
  ps_obj@sam_data$sample_type = ps_obj@sam_data$sample_type %>% tolower()
  ps_obj@tax_table@.Data = cbind(ps_obj@tax_table@.Data, taxa_names(ps_obj)); colnames(ps_obj@tax_table@.Data)[ncol(ps_obj@tax_table@.Data)] = "asv"
  ps_blank_obj = ps_obj %>% subset_samples(sample_type != "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  ASVs_in_Blanks = taxa_names(ps_blank_obj)

  if(method == "max_v"){
    Extraction_neg_max_vec <- apply(ps_blank_obj %>% pstoveg_otu %>% t() %>% as.data.frame, 1, max) %>% as.vector()
    names(Extraction_neg_max_vec) = taxa_names(ps_blank_obj)
    Extractiondf = ps_obj %>% subset_taxa(asv %in% ASVs_in_Blanks) %>% pstoveg_otu %>% as.data.frame()
    Extractiondf = sweep(Extractiondf,MARGIN=2,Extraction_neg_max_vec,FUN="-")
    Extractiondf <- replace(Extractiondf, Extractiondf < 0, 0)
    new_df = ps_obj %>% subset_taxa(asv %ni% ASVs_in_Blanks) %>% pstoveg_otu()
    new_df = cbind(new_df, Extractiondf)
    ps_trimmed_obj = ps_obj
    ps_trimmed_obj@otu_table = otu_table(new_df, taxa_are_rows=FALSE)
    ps_trimmed_obj = ps_trimmed_obj %>% subset_samples(sample_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  } else{
    ps_trimmed_obj = ps_obj %>% subset_samples(sample_type == "sample") %>% subset_taxa(asv %ni% ASVs_in_Blanks)
  }
  # Printing results
  ntaxa_before = ps_obj %>% subset_samples(sample_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% ntaxa()
  ntaxa_after = ntaxa(ps_trimmed_obj)

  cat("\n")
  cat(paste0("Contamination removal outcome using ",method),"\n")
  cat(paste0("Number of ASVs totally removed: ",ntaxa_before - ntaxa_after))
  cat("\n")
  cat(paste0("Percent of ASVs removed: ",round((1 - (ntaxa_after / ntaxa_before)) * 100,2), " %"))
  cat("\n")

  reads_before = ps_obj %>% subset_samples(sample_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% sample_sums() %>% sum()
  reads_after = ps_trimmed_obj %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% sample_sums() %>% sum()

  cat(paste0("Total number of reads  removed: ",reads_before - reads_after))
  cat("\nPercent of reads removed: ",paste0(round((1 - (reads_after / reads_before)) * 100,2)," %"))
  cat("\n")
  return(ps_trimmed_obj)
}
