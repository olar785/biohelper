#' pstoveg_sample
#'
#' @description
#' This function from the theseus R package takes in a phyloseq object and returns a metadata table compatible with the vegan R package.
#' If you use this function for a publication, please cite the creators of the theseus R package (citation('theseus')).
#
#'
#' @param
#' PS             (required) a phyloseq object
#' @export
#' @examples
#' pstoveg_sample(PS = ps_test_data)

pstoveg_sample = function (PS)
{
  SAMP <- phyloseq::sample_data(PS)
  return(as(SAMP, "data.frame"))
}
