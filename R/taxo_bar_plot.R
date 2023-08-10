#' taxo_bar_plot
#'
#' @description
#' This function takes in a phyloseq object with taxonomy and returns a bar plot with colours and shades per
#' taxa. The output is not directly customizable but can be re-used to add facets afterward.
#' Sample can also be ordered using the s_order argument.
#'
#' @param
#' rank1                Taxonomic rank to be associated with colours (e.g. Phylum)
#' @param
#' rank2                Taxonomic rank to be associated with shades of colours (e.g. Family)
#' @param
#' colors               Colors to be used
#' @param
#' n_rank2              Number of shades to use per colour. Taxa which rank below that number will be assigned to "Other".
#' @param
#' alpha_num_ordering   Whether to order samples alpha-numerically (default == F)
#'
#' @export
#' @examples
#' colors = c("cyan", "palegreen", "yellow", "deeppink ", "white", "dodgerblue", "lightsalmon")\cr
#' ps_test_data_t = ps_test_data %>% tax_glom('Family') %>% microbiome::transform(transform = "compositional")\cr
#' p1 = taxo_bar_plot(ps_test_data_t, rank1 = "Phylum", rank2 = "Family", colors = colors,  f = "extraction_method")\cr
#' p1 + facet_wrap(extraction_method~., drop = TRUE, scale="free", nrow = 1) + ggtitle("Taxonomic composition per extraction method")

taxo_bar_plot = function(ps_obj, rank1 = "Phylum", rank2 = "Family", n_rank1 = NA, n_rank2 = 6, colors = c("cyan", "palegreen", "yellow", "deeppink ", "white", "dodgerblue", "lightsalmon"), label = NA, alpha_num_ordering=F){
  if(is.na(label)){
    label = "Sample"
  }

  n_rank1 = if(is.na(n_rank1)){
    n_rank1 = length(colors)-1
  }
  # Creating dataframe
  dfn = speedyseq::psmelt(ps_obj) %>% base::replace(is.na(.), "Others")
  # Removing taxa with zero abundance
  dff = dfn %>% filter(Abundance >0)
  # Converting all factors to character strings
  dff = dff %>% mutate_if(is.factor, as.character)
  # Changing rare rank1 to 'Others'
  ranks = colnames(ps_obj@tax_table@.Data %>% remove_empty("cols"))
  #dff[,rank1][grepl("_X", dff[,rank1], ignore.case=FALSE)] <- "Others"
  temp = dff %>% dplyr::group_by_at(which(colnames(dff)==ranks[1]):which(colnames(dff)==rank1)) %>% dplyr::summarise(rank1_Sum = sum(Abundance))
  temp_2 = temp %>% arrange(desc(rank1_Sum))
  temp_2 = temp_2 %>% head(n_rank1)
  list_of_rare_rank1 = temp[(temp[,rank1] %>% as_vector() %>% unname()) %ni% (temp_2[,rank1] %>% as_vector() %>% unname()),rank1] %>% as_vector() %>% unname()
  dff[dff[,rank1] %in% list_of_rare_rank1,rank1] = "Others"
  # Changing rare families to 'Others' including those with _X...
  dff[,rank2][grepl("_X", dff[,rank2], ignore.case=FALSE)] <- "Others"
  temp = dff %>% dplyr::group_by_at(c(which(colnames(dff)=="OTU"), which(colnames(dff)==ranks[1]):which(colnames(dff)==ranks[length(ranks)]))) %>% dplyr::summarise(ran2_Sum = sum(Abundance)) %>% as.data.frame()
  temp_3 = temp %>%
    dplyr::group_by_at(c(which(colnames(temp) == ranks[1]):which(colnames(temp) == rank1))) %>%
    dplyr::arrange(desc(ran2_Sum), .by_group = T) %>%
    as.data.frame() %>%
    slice_head(n = 5,by = rank1)
  list_of_rare_rank2 = temp[(temp[,rank2] %>% as_vector() %>% unname()) %ni% (temp_3[,rank2] %>% as_vector() %>% unname()),]$OTU
  list_of_rare_rank2 = c(list_of_rare_rank2,temp[temp[,which(colnames(temp)==rank2)]=="NA",]$OTU)
  dff[dff$OTU %in% list_of_rare_rank2,rank2] = "Others"
  #darken each color n times in increments of steps towards black
  ExpandColors <- function(colors, n, steps = 11){
    if(n <= steps){
      suppressWarnings({
        sapply(colors, function(x){colorRampPalette(c(x, "#000000"))(steps)}) %>%
          as.data.frame() %>%
          filter(row_number() <= n) %>%
          gather(key = original.color, value = expanded.color)
      })
    }else{
      warning("Select n < steps!")
    }
  }
  r1 = unique(dff[,rank1]) %>% as.character() %>% sort()
  names(r1) = r1
  color_list = list()
  i = 1
  for (j in seq(from = 1, to = length(colors))) {
    getPalette = colorRampPalette(colors[i], alpha = T)
    taxaList = sort(unique(dff[(dff[,rank1] %>% as.vector() %>% unname()) == names(r1[i]),rank2] %>% as.vector() %>% unname()))
    color_list[i] = list(ExpandColors(colors[i], n = length(taxaList)))
    r1[i] = list(as.data.frame(color_list[i])$expanded.color)
    if (!is.null(r1[[i]])) {
      names(r1[[i]]) = taxaList
    }
    i = i+1
  }
  r1Palette = unlist(r1)
  # Adding r1 name to r2 name
  dff[,rank2] = paste0(dff[,rank1] %>% as.vector() %>% unname(),".",dff[,rank2] %>% as.vector() %>% unname())
  dff = dff[order(dff[,rank2]),]
  labls = unique(dff[,rank2])
  names(labls) = seq(length(labls))
  # Adding a labbel column
  func = function(rank2){
    return(names(labls[labls == rank2]))
  }
  dff$labls = sapply(dff[,rank2], func)
  dft = dff %>% dplyr::group_by_at(which(colnames(dff) %in% c("Sample",rank2))) %>% dplyr::summarise(Abund=sum(Abundance)*100)
  dff = dff[match(dft$Sample, dff$Sample),]
  dft = cbind(dft, dff[,colnames(dff)%ni%colnames(dft)])
  if(alpha_num_ordering==T){
    dft$Sample = factor(dft$Sample, levels = gtools::mixedsort(dft$Sample %>% unique()))
  }

  # Making plot
  return(ggplot(data = dft, aes(Sample, Abund, fill = pull(dft, rank2), label = pull(dft, label))) +
           geom_bar(stat="identity",color="black") +
           theme_bw() +
           ylab("Relative abundance (%)") +
           scale_fill_manual(values=r1Palette) +
           theme_classic() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
           guides(fill=guide_legend(title=paste0(str_to_title(rank1),"; ",str_to_title(rank2))))
  )
}
