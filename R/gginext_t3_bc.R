#' gginext_t3_bc
#'
#' @description
#' This function takes in an iNEXT object and create a ggiNEXT plot of type 3 (Diversity per sample coverage) that contains a vertical line at base coverage (defined as the highest coverage value between minimum extrapolated values and maximum interpolated values)
#' @param
#' inext          iNEXT onject
#'
#' @export
#' @examples
#' inext_test = inext_ifreq_wrapper(pst = ps_test_data,  grp = "biome", q = 0, knots = 100, nboot = 200)
#' gginext_t3_bc(inext = inext_test)

gginext_t3_bc = function(inext){
  temp = inext$iNextEst$coverage_based %>% dplyr::filter(Method=="Extrapolation") %>% dplyr::group_by(Assemblage) %>% top_n(1, SC)
  Base_coverage = min(temp$SC)
  iNextp = iNEXT::ggiNEXT(inext, type=3) + geom_vline(xintercept= Base_coverage, linetype="dashed", color = "black")
  return(iNextp)
}
