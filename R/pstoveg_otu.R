#' pstoveg_otu
#'
#' @description
#' This function from the theseus R package takes in a phyloseq object and returns a OTU/ASV table compatible with the vegan R package.
#' If you use this function for a publication, please cite the creators of the theseus R package (citation('theseus')).
#'
#' @param
#' PS             (required) a phyloseq object
#' @export
#' @examples
#' pstoveg_otu(PS = ps_test_data)


pstoveg_otu = function (PS)
{
  OTU <- phyloseq::otu_table(PS)
  if (phyloseq::taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
