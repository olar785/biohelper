#' Load a taxo_normalisation
#'
#' This function performs normalisation of taxonomic assignments of 18S and COI data
#' identified via silva, midori or any other reference databases, using NCBI curated
#' classification and nomenclatureand database. It can also provide further taxonomic
#' information of additional desired ranks.
#'
#' The first column must contain ASV/OTU identifiers as the example below:
#' ASV     Kingdom   Phylum   etc
#' ASV_1   Eukaryota   Arthropoda
#' ASV_2   Eukaryota   Chordata
#' ASV_3   Eukaryota   Chordata
#'
#' Importantly, this function requires the download of NCBI taxonomic database (~3 GB).
#' This takes few minutes to install using the following command:
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package

#'
#' @export
#' @examples
#' taxo_normalisation(taxo_dataframe, sqlFile = sqlFile, desired_ranks = c("Family", "Genus", "Species"))
taxo_normalisation = function(df, sqlFile, ranks){
  df[df==""]<-NA
  paternsToRemove = c("^.+_environmental.+|environmental_.+|uncultured_.+|_sp\\..+|_sp.|_sp.+| sp\\..+| sp.| sp.+|_\\(.+|^.+_metagenome|_cf.")
  df = df %>% mutate_all(list(~str_replace(.,paternsToRemove, ""))) %>% mutate_all(list(~na_if(.,"")))
  if("species" %in% tolower(colnames(df))){
    df$Species = paste0(sapply(strsplit(df$Species,"_"), `[`, 1)," ", sapply(strsplit(df$Species,"_"), `[`, 2))
  }
  df = df %>% mutate_all(list(~str_replace(.,"NA NA| NA", ""))) %>% mutate(across(everything(), gsub, pattern = "_", replacement = " ")) %>% mutate_all(list(~na_if(.,"")))
  ranks_indexes = which(tolower(colnames(df)) %in% ranks)
  non_taxo_ranks = c("otu","otus","OTU","OTUs","OTUS","asv","asvs","ASVs","ASV","ASVS","nR")
  rpt_indexes = max.col(!is.na(df[colnames(df)%ni%non_taxo_ranks]), "last")
  taxa = unlist(lapply(1:length(rpt_indexes), function(x) df[x, rpt_indexes[x]]))
  res_df = data.frame("ASV" = rownames(df), "rpt_indexes" = rpt_indexes, "taxa" = taxa)
  res_df$id = getId(taxa = res_df$taxa, sqlFile = sqlFile, onlyScientific = TRUE)
  length(res_df[which(is.na(res_df$id)),]$id)
  r = 1
  # Deals with NA ids
  while (r<length(ranks_indexes) & any(is.na(res_df$id))) {
    df_temp = df[which(is.na(res_df$id)),]
    res_df_temp = res_df[which(is.na(res_df$id)),]
    rpt_indexes = max.col(!is.na(df_temp[colnames(df_temp)%ni%non_taxo_ranks]), "last") - r
    rpt_indexes = pmax(rpt_indexes,1) # makes sure to have no negative or 0 values
    taxa = unlist(lapply(1:length(rpt_indexes), function(x) df_temp[x, rpt_indexes[x]]))
    res_df_temp = data.frame("ASV" = rownames(df_temp), "rpt_indexes" = rpt_indexes, "taxa" = taxa)
    res_df[which(is.na(res_df$id)),]$id = getId(taxa = res_df_temp$taxa, sqlFile = sqlFile, onlyScientific = TRUE)
    length(res_df[which(is.na(res_df$id)),]$id)
    r = r + 1
  }
  length(res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id)
  # Deals with multiple ids
  n = 1
  number_ids = max(str_count(res_df$id, pattern = ","),na.rm = T) + 1
  while (n<=number_ids & any(str_detect(res_df$id, ",", negate = FALSE))) {
    df_temp = df[which(str_detect(res_df$id, ",", negate = FALSE)),]
    res_df_temp = res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]

    rpt_indexes = max.col(!is.na(df_temp[colnames(df_temp)%ni%non_taxo_ranks]), "last") - 1
    rpt_indexes = pmax(rpt_indexes,1) # makes sure to have no negative or 0 values
    res_df_temp$p_taxa = unlist(lapply(1:length(rpt_indexes), function(x) df_temp[x, rpt_indexes[x]]))

    taxa = unlist(lapply(1:length(rpt_indexes), function(x) df_temp[x, rpt_indexes[x]]))
    test = getTaxonomy(res_df_temp$id %>% strsplit( "," ) %>% sapply( "[", n ), sqlFile, desiredTaxa = ranks) %>% as.data.frame()
    for (i in 1:nrow(res_df_temp)) {
      res_df_temp$id[i] = ifelse(res_df_temp$p_taxa[i] %in% test[i,], res_df_temp$id[i] %>% strsplit( "," ) %>% sapply( "[", n ),res_df_temp$id[i])
    }
    res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id = res_df_temp$id
    length(res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id)
    n = n + 1
  }
  length(res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id)
  # Deals with multiple ids again but at higher level
  n = 1
  number_ids = max(str_count(res_df$id, pattern = ","),na.rm = T) + 1
  if(!is.na(any(str_detect(res_df$id, ",", negate = FALSE))))
    while (n<=number_ids & any(str_detect(res_df$id, ",", negate = FALSE))) {
      df_temp = df[which(str_detect(res_df$id, ",", negate = FALSE)),]
      res_df_temp = res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]

      rpt_indexes = max.col(!is.na(df_temp[colnames(df_temp)%ni%non_taxo_ranks]), "last") - 2
      rpt_indexes = pmax(rpt_indexes,1) # makes sure to have no negative or 0 values
      res_df_temp$p_taxa = unlist(lapply(1:length(rpt_indexes), function(x) df_temp[x, rpt_indexes[x]]))

      taxa = unlist(lapply(1:length(rpt_indexes), function(x) df_temp[x, rpt_indexes[x]]))
      test = getTaxonomy(res_df_temp$id %>% strsplit( "," ) %>% sapply( "[", n ), sqlFile, desiredTaxa = ranks) %>% as.data.frame()
      for (i in 1:nrow(res_df_temp)) {
        res_df_temp$id[i] = ifelse(res_df_temp$p_taxa[i] %in% test[i,], res_df_temp$id[i] %>% strsplit( "," ) %>% sapply( "[", n ),res_df_temp$id[i])
      }
      res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id = res_df_temp$id
      length(res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id)
      n = n + 1
    }
  #print(length(res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id))
  res_df[ranks] = getTaxonomy(res_df$id, sqlFile, desiredTaxa = ranks)
  res_df = res_df[,colnames(res_df) %in% c(ranks,non_taxo_ranks)]
  res_df$superkingdom = res_df$superkingdom %>% replace_na("Unknown")
  res_df = res_df %>% column_to_rownames("ASV")
  res_df=cbind(res_df,df[,colnames(df) %in% non_taxo_ranks,drop = FALSE])
  return(res_df)
}
