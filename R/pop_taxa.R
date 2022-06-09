#' Loads pop_taxa
#'
#' This function takes in a phyloseq object and a list of taxa to be removed.
#' Credit to joey711 (https://github.com/joey711/phyloseq/issues/652)
#'
#' ps_obj = Phyloseq object to trim.
#' undesiredTaxa = names of taxa to be removed (generally ASV or OTU id), as long as it is in taxa_names() of the phyloseq object.
#'
#' @export
#' @examples
#' undesiredTaxa = c("3b8f1e9447e2b3f55113dcd5d04eb152", "3efe4018e327153524ce8feb1db016ea", "3efe4018e327153524ce8feb1db016ea")
#' ps_test_data = ps_test_data %>% pop_taxa(undesiredTaxa)

pop_taxa = function(ps_obj, undesiredTaxa){
  allTaxa = taxa_names(ps_obj)
  myTaxa <- allTaxa[!(allTaxa %in% undesiredTaxa)]
  return(prune_taxa(myTaxa, ps_obj))
}
