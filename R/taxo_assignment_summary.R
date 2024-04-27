#' taxo_assignment_summary
#'
#' This function takes in a phyloseq object with taxonomy or a dataframe of taxonomic assignment and returns the percent of assigned taxa per taxonomic rank, and the mean number and standard deviation of taxonomic rank assigned.
#'
#' @param
#' obj A phyloseq object with taxonomy or a dataframe of taxonomic assignment
#' @param
#' ranks  Taxonomic ranks of interest from the taxonomic assignment. By default, the function uses colnames of the tax_table of the phyloseq object or of the dataframe, depending on what is provided.
#' print_only If TRUE (default), results are printed only. If FALSE, a list containing 1) a dataframe of percent of assigned taxa per taxonomic rank, and 2) a dataframe of assignment per taxonomic rank per taxa (0 = unassigned, 1 = assigned) are returned.
#'
#' @export
#' @examples
#' taxo_assignment_summary(obj = ps_test_data, ranks = c("Kingdom","Phylum","Class","Order","Family","Genus"))
#'
#' # Testing if taxonomic assignation at a specific rank is different between two groups:
#' ps_manual_ta = ps_test_data %>% subset_samples(extraction_method == "manual") %>% phyloseq::filter_taxa(function(x) sum(x) >0, TRUE) %>% taxo_assignment_summary(print_only = F)
#' ps_robot_ta = ps_test_data %>% subset_samples(extraction_method == "robot") %>% phyloseq::filter_taxa(function(x) sum(x) >0, TRUE) %>% taxo_assignment_summary(print_only = F)
#' t.test(x = ps_manual_ta$df_assignment_per_taxaNrank %>% dplyr::filter(taxonomic_rank == "Family") %>% dplyr::pull(assignment),
#'        y = ps_robot_ta$df_assignment_per_taxaNrank %>% dplyr::filter(taxonomic_rank == "Family") %>% dplyr::pull(assignment))

taxo_assignment_summary = function(obj, ranks = NA, print_only = T){
  if(any(class(obj) %in% "phyloseq")){
    df = obj@tax_table@.Data %>% as.data.table()
  }else{
    df = obj %>% as.data.table()
  }
  if(any(is.na(ranks))){
    ranks = colnames(df)[tolower(colnames(df)) %ni% c("asv", "otu")]
  }
  df[df==""]<-NA
  df$nRb = length(ranks) - rowSums(is.na(df %>% dplyr::select(all_of(ranks))))
  summary_overall = df %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  temp_summary_per_rank = df  %>%
    pivot_longer(cols = ranks, names_to = "taxonomic_rank", values_to = "taxa") %>%
    dplyr::mutate(assignment = ifelse(is.na(taxa), 0, 1)) %>%
    dplyr::select(-nRb)
  summary_per_rank = temp_summary_per_rank %>%
    dplyr::group_by(taxonomic_rank) %>%
    dplyr::summarise(percent_assigned = round(sum(assignment) / nrow(df),2)*100) %>%
    dplyr::arrange(desc(percent_assigned))

  if(print_only != T){
    lt <- sapply(c("df_summary_per_rank","df_assignment_per_taxaNrank"),function(x) NULL)
    lt[["df_summary_per_rank"]] = summary_per_rank
    lt[["df_assignment_per_taxaNrank"]] = temp_summary_per_rank
    return(lt)
  }else{
    print(summary_per_rank)
  }
  cat("\nMean number of assigned taxonomic rank: ",summary_overall$mean %>% as.numeric(),"( out of",length(ranks),")","\nStandard deviation: ",summary_overall$sd %>% as.numeric(),"\n\n")
}
