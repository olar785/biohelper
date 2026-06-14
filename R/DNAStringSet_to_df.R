#' DNAStringSet_to_df
#'
#' This function takes in a DNAStringSet object and return a dataframe.
#'
#' @param dss A DNAStringSet object.
#'
#' @export
#' @examples
#' DNAStringSet_to_df(ps_test_data@refseq)
#'
DNAStringSet_to_df = function(dss){
  return(data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss)))
}
