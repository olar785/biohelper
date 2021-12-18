#' Loads DNAStringSet_to_df
#'
#' This function takes in a DNAStringSet object and return a dataframe.
#'
#' @export
DNAStringSet_to_df = function(dss){
  return(data.frame(width=width(dss), seq=as.character(dss), names=names(dss)))
}
