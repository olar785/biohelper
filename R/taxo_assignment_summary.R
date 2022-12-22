#' taxo_assingment_summary
#'
#' This function takes in a dataframe of taxonomic assignment and returns the mean number and standard deviation of taxonomic ranks assigned.
#'
#' @param
#' df Dataframe of taxonomic assignment
#' @param
#' ranks  Taxonomic ranks in the dataframe
#'
#' @export
#' @examples
#' df = ps_test_data@tax_table %>% as.data.frame()
#' taxo_assingment_summary(df = df, ranks = c("Kingdom","Phylum","CLass","Order","Family","Genus"))

taxo_assignment_summary = function(df, ranks){
  df = df %>%
    as.data.table() %>%
    na_if('')
  df$nRb = length(ranks) - rowSums(is.na(df %>% dplyr::select(ranks)))
  temp_summary = df %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric(),"\n\n")
}
