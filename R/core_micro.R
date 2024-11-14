#' core_micro
#'
#' @description
#' This function takes in a phyloseq object and subset core taxa per a specified group as opposed to subsetting by the overall data.
#' Two methods  currently available are 1) the 'core' function of the microbiome R package (Lahti et al., 2019) and 2) subsetting by Cumulative Relative ABundance (crab).
#' The latter is defined by setting a divider to the maximum cumulative abundance, which will be used as a minimum threshold value to keep taxa.
#' If you use this function for a publication, please cite the creators of the microbiome R package (citation('microbiome')).
#'
#' @param
#' pst            Phyloseq object
#' @param
#' grp            Group to perform the core microbiome subsetting by.
#' @param
#' method         Subsetting methodology. Two methods are currently available: 'core' and 'crab'
#' @param
#' n          	  Divider value for the 'crab' method. To keep ASVs with cumulative relative abundance within two orders of magnitude of the maximum cumulative relative abundance for example, n would be set at 100.
#' @param
#' d        	  Detection threshold for the 'core' function.
#' @param
#' p        	  Prevalence threshold for the 'core' function.
#' @param
#' ...        	  Other arguments to pass to the 'core' function.
#'
#' @export
#' @examples
#' core_micro(pst = ps_test_data,  grp = "biome", method = "crab", n = 10)


core_micro = function(pst, grp, method = NULL, n = NULL, d = NULL, p = NULL, ...){
  pst@phy_tree=NULL
  mt = pst %>% pstoveg_sample()
  mylist.names = mt %>% dplyr::pull(grp) %>% unique()
  pstmerged = pst %>%
    microbiome::transform("compositional") %>%
    merge_samples(group = grp, fun = sum)
  pstmerged@sam_data$sid = sample_names(pstmerged)

  mylist_taxa <- sapply(mylist.names,function(x) NULL)
  mylist_pst <- sapply(mylist.names,function(x) NULL)
  pstmergedmelt = biohelper::psmelt(pstmerged)
  pstmelt = biohelper::psmelt(pst)
  taxn = pst@tax_table@.Data %>% colnames()
  seqs = pst@refseq

  if(method == "abun")
    for (i in names(mylist_taxa)) {
      mylist_taxa[[i]] = pstmergedmelt %>% dplyr::filter(Sample == i)
      s = mylist_taxa[[i]]$Abundance %>% sum() / n
      mylist_taxa[[i]] = mylist_taxa[[i]] %>% dplyr::filter(Abundance >= s) %>% dplyr::pull(OTU)
      mylist_pst[[i]] = pstmelt %>% dplyr::filter(!!as.symbol(grp) ==i & OTU %in% mylist_taxa[[i]])
      # Recreating phyloseq object
      otu_tab = spread(mylist_pst[[i]] %>% dplyr::select(c(OTU,Sample,Abundance)), key = Sample, value = Abundance) %>%
        column_to_rownames("OTU") %>%
        otu_table(taxa_are_rows = T)
      tax_tab = mylist_pst[[i]] %>%
        dplyr::select(c(OTU,taxn)) %>%
        distinct() %>%
        column_to_rownames("OTU") %>%
        as.matrix() %>%
        tax_table()
      sam_tab = mt %>% dplyr::filter(!!as.symbol(grp) == i) %>% sample_data()
      seqs_tab = seqs[seqs@ranges@NAMES %in% rownames(tax_tab)] %>% refseq()
      mylist_pst[[i]] = merge_phyloseq(otu_tab, tax_tab, sam_tab, seqs_tab)
    }

  if(method == "core"){
    for (i in names(mylist_taxa)) {
      mylist_pst[[i]] = pstmelt %>% dplyr::filter(!!as.symbol(grp) == i)
      # Recreating phyloseq object
      otu_tab = spread(mylist_pst[[i]] %>% dplyr::select(c(OTU,Sample,Abundance)), key = Sample, value = Abundance) %>%
        column_to_rownames("OTU") %>%
        otu_table(taxa_are_rows = T)
      tax_tab = mylist_pst[[i]] %>%
        dplyr::select(c(OTU,taxn)) %>%
        distinct() %>%
        column_to_rownames("OTU") %>%
        as.matrix() %>%
        tax_table()
      sam_tab = mt %>% dplyr::filter(!!as.symbol(grp) == i) %>% sample_data()
      seqs_tab = seqs[seqs@ranges@NAMES %in% rownames(tax_tab)] %>% refseq()
      mylist_pst[[i]] = merge_phyloseq(otu_tab, tax_tab, sam_tab, seqs_tab)
      # Filtering with microbiome package
      mylist_pst[[i]] <- microbiome::core(mylist_pst[[i]], detection = d, prevalence = p, ...)
    }}
  mergers <- list()
  for (sam in mylist_pst) {
    mergers = merge_phyloseq(mergers,sam)
  }
  return(mergers)
}
