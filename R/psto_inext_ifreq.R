#' psto_inext_ifreq
#'
#' @description
#' This function takes in a phyloseq object and perform iNEXT from the iNEXT R package (Hsieh et al., 2022) on incidence frequency data (Chao2 method).
#' Briefly, the data is split into a list of dataframes based on the grouping variable and frequency of incidence per group is calculated.
#' If you use this function for a publication, please cite the creators of the iNEXT R package (citation('iNEXT')).
#'
#' @param
#' pst            Phyloseq object
#' @param
#' q              A number or vector specifying the diversity order(s) of Hill numbers.
#' @param
#' grp            Group to compare
#' @param
#' ...        	  Other arguments to pass to the 'iNEXT' function.
#'
#' @export
#' @examples
#' inext_ifreq_wrapper(pst = ps_test_data,  grp = "biome", q = 0, knots = 100, nboot = 200)

psto_inext_ifreq = function (pst, q = 0, grp, ...)
{
  mt = pst %>% pstoveg_sample()
  mylist.names = mt %>% dplyr::pull(grp) %>% unique()
  mylist_pst <- sapply(mylist.names, function(x) NULL)
  pstmelt = speedyseq::psmelt(pst)
  for (i in mylist.names) {
    mylist_pst[[i]] = pstmelt %>% dplyr::filter(!!as.symbol(grp) == i)
    otu_tab = spread(mylist_pst[[i]] %>% dplyr::select(c(OTU, Sample, Abundance)), key = Sample, value = Abundance) %>%
      tibble::column_to_rownames("OTU") %>% phyloseq::otu_table(taxa_are_rows = T)
    sam_tab = mt %>% dplyr::filter(!!as.symbol(grp) == i) %>%
      phyloseq::sample_data()
    mylist_pst[[i]] = phyloseq::merge_phyloseq(otu_tab, sam_tab) %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
    mylist_pst[[i]] = mylist_pst[[i]] %>% microbiome::transform(transform = "pa") %>%
      pstoveg_otu() %>% t() %>% as.data.frame()
  }
  B = lapply(mylist_pst, iNEXT::as.incfreq)
  inext_q0 <- iNEXT::iNEXT(B, q = q, datatype = "incidence_freq", ...)
  return(inext_q0)
}




