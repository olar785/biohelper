#' Loads NAs_to_Xs
#'
#' This function replaces NA entries in a taxonomic table by the previous (higher)
#' taxonomic information and adds Xs according to the number of NA entries for each organism.
#'
#' The first column must contain ASV/OTU identifiers as the example below:
#' ASV     Kingdom     Phylum   Genus
#' ASV_1   Eukaryota
#' ASV_2   Eukaryota   Chordata  Mus
#' ASV_3   Eukaryota   Chordata
#'
#' @export
#' @examples
#' NAs_to_Xs(ps_test_data@tax_table@.Data %>% as.data.frame())

NAs_to_Xs = function(df){
  df[df==""]<-NA
  non_taxo_ranks = c("otu","otus","OTU","OTUs","OTUS","asv","asvs","ASVs","ASV","ASVS")
  ranks_indexes = which(tolower(colnames(df)) %ni% non_taxo_ranks)
  df[is.na(df[,min(ranks_indexes)]) | df[,min(ranks_indexes)]=="" ,min(ranks_indexes)] = "Unknown"
  for (i in (min(ranks_indexes)+1):max(ranks_indexes)) {
    df[which(is.na(df[,i])),i] = paste0(df[which(is.na(df[,i])),(i-1)],"_X")
  }
  df = df %>% mutate_all(list(~str_replace(., "_", "")))
  df = data.frame(lapply(df, function(x) {gsub("_", "", x)}))
  df = data.frame(lapply(df, function(x) {gsub("(^.)+(X+)", "\\1_\\2", x)}))
  df = data.frame(lapply(df, function(x) {gsub("(X+$)", "_\\1", x)}))
  return(df)
}
