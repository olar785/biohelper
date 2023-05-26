#' upset_wrapper
#'
#' @description
#' This function is a wrapper to create an 'upset' figure from a phyloseq object and using the UpSetR package by Gehlenborg (2019).
#' It returns a dataframe for UpsetR, the plot obtained directly from the 'upset' function and a ggplot2 object obtained from the cowplot::plot_grid function for conveniency.
#'
#' @param
#' pst            Phyloseq object
#' @param
#' grp            Group to compare shared and unique taxon.
#' @param
#' order.by       How the intersections in the matrix should be ordered by. Options include frequency (entered as "freq"), degree, or both in any order.
#' @param
#' sets.x.label   The x-axis label of the set size bar plot.
#' @param
#' rel_widths     Optional argument for the plot_grid function: Numerical vector of relative columns widths. For example, in a two-column grid, rel_widths = c(2, 1) would make the first column twice as wide as the second column.
#' @param
#' rel_heights    Optional argument for the plot_grid function: Numerical vector of relative rows heights. Works just as rel_widths does, but for rows rather than columns.
#' @param
#' ...        	  Other arguments to pass to the 'upset' function.
#'
#' @export
#' @examples
#' upset_wrapper(pst = ps_test_data, grp = "biome")

upset_wrapper = function(pst, grp, order.by = "freq", sets.x.label = "Taxa richness", nrow = NULL, align = NULL, rel_heights = NULL, rel_widths = NULL,...){
  mylist <- sapply(c("upset_df","upset_plot","upset_plotv2"),function(x) NULL)
  upset = pst %>% 
    phyloseq::filter_taxa(function(x) sum(x) >0, TRUE) %>%
    merge_samples(grp) %>% 
    microbiome::transform("pa")
  mylist[[1]] = upset %>% pstoveg_otu() %>% t() %>% as.data.frame()
  ups = upset(mylist[[1]], order.by = order.by, sets.x.label = sets.x.label, nsets = nrow(mylist[[1]]), ...)
  mylist[[3]] = cowplot::plot_grid(NULL, mylist[[2]]$Main_bar, mylist[[2]]$Sizes, mylist[[2]]$Matrix,nrow= 2, align='hv', rel_heights = rel_heights, rel_widths = rel_widths)
  return(mylist)
}