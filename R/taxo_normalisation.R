#' Loads taxo_normalisation
#'
#' @description
#' This function performs normalisation of taxonomic assignments of 18S and COI data
#' identified via silva, midori or any other reference databases, using NCBI curated
#' classification and nomenclatureand database. It can also provide further taxonomic
#' information of additional desired ranks. It can either take in a taxonomy table
#' (e.g. from Qiime2) or a phyloseq object.
#'
#' Importantly, this function requires the download of NCBI taxonomic database (~65 GB).
#'
#' This takes few minutes to install using the following command:\cr
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package
#'
#' @param
#' obj = Either a dataframe containing taxonomic information (e.g. output of Qiime2 or the assignTaxonomy from dada2) or a phyloseq object
#' @param
#' sqlFile = Path to the local NCBI taxonomy db
#' @param
#' ranks = Ranks to return
#' @param
#' keepSAR = Keep the SAR assignment from the input data, which is not a valid group in NCBI taxonomy db (default is FALSE)
#'
#' @export
#' @examples
#' data("ps_test_data")
#' taxo_normalisation(ps_test_data, sqlFile = 'accessionTaxa.sql', ranks = c("Superkingdom", "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus"))
#'

taxo_normalisation = function(obj, sqlFile, ranks, keepSAR = F){

  if("phyloseq" %in% class(obj)){
    df = obj@tax_table@.Data %>% as.data.frame()
  } else{
    df = obj %>% as.data.frame()
  }

  ranks = tolower(ranks)
  colnames(df) = tolower(colnames(df))
  taxa_id = c("otu","otus","asv","asvs","feature_id","feature.id","nR")

  if("taxon" %in% colnames(df)){
    df$taxon = gsub("D_\\d+__","",taxo$taxon)
    df = df %>% dplyr::select(c(colnames(df)[which(colnames(df)%in%c(taxa_id,"taxon"))])) %>% splitstackshape::cSplit('taxon', sep = ';')
  }else if("taxonomy" %in% colnames(df)){
    df = df %>% dplyr::select(c(colnames(df)[which(colnames(df)%in%c(taxa_id,"taxonomy"))])) %>% splitstackshape::cSplit('taxonomy', sep = ';')
  }
  df = df %>% mutate_if(is.factor, as.character)
  df[is.na(df)]<-""

  if(any(taxa_id %in% colnames(df))){
    taxa_id_df = which(colnames(df) %in% taxa_id)%>% as.numeric()
    colnames(df)[taxa_id_df] = "feature_id"
    if(any(duplicated(df))){
      cat("Duplicated rows! Keeping the first instance only.\n")
      df = df[!duplicated(df), ]
      rownames( df ) <- NULL
      df = df %>% column_to_rownames("feature_id")
    }else{
      rownames( df ) <- NULL
      df = df %>% column_to_rownames("feature_id")
    }
  }

  paternsToRemove = c("^.+_environmental.+|environmental_.+|uncultured_.+|_sp\\..+|_sp.|_sp.+| sp\\..+| sp.| sp.+|_\\(.+|^.+_metagenome|_cf.")
  df = df %>% mutate_all(list(~str_replace(.,paternsToRemove, ""))) %>% mutate_all(list(~na_if(.,"")))
  if("species" %in% colnames(df)){
    df$Species = paste0(sapply(strsplit(df$species,"_"), `[`, 1)," ", sapply(strsplit(df$species,"_"), `[`, 2))
  }

  df = df %>% mutate_all(list(~str_replace(.,"NA NA| NA", ""))) %>% mutate(across(everything(), gsub, pattern = "_", replacement = " ")) %>% mutate_all(list(~na_if(.,"")))
  ranks_indexes = which(colnames(df) %ni% c("otu","otus","asv","asvs","feature_id","feature.id","nR"))
  non_taxo_ranks = c("otu","otus","asv","asvs","feature_id","feature.id","nR")
  rpt_indexes = max.col(!is.na(df[ranks_indexes]), "last")

  taxa = unlist(lapply(1:length(rpt_indexes), function(x) df[x, rpt_indexes[x]]))
  res_df = data.frame("feature_id" = rownames(df), "rpt_indexes" = rpt_indexes, "taxa" = taxa)
  res_df$id = getId(taxa = res_df$taxa, sqlFile = sqlFile, onlyScientific = TRUE)
  r = 1
  # Deals with NA ids
  while (r<length(ranks_indexes) & any(is.na(res_df$id))) {
    df_temp = df[which(is.na(res_df$id)),]
    res_df_temp = res_df[which(is.na(res_df$id)),]
    rpt_indexes = max.col(!is.na(df_temp[colnames(df_temp)%ni%non_taxo_ranks]), "last") - r
    rpt_indexes = pmax(rpt_indexes,1) # makes sure to have no negative or 0 values
    taxa = unlist(lapply(1:length(rpt_indexes), function(x) df_temp[x, rpt_indexes[x]]))
    res_df_temp = data.frame("feature_id" = rownames(df_temp), "rpt_indexes" = rpt_indexes, "taxa" = taxa)
    id = getId(taxa = res_df_temp$taxa, sqlFile = sqlFile, onlyScientific = TRUE)
    res_df[which(is.na(res_df$id)),]$id = ifelse(!is.na(id),is,"")
    r = r + 1
  }
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
    number_ids = max(str_count(res_df$id, pattern = ","),na.rm = T) #####
  }
  #length(res_df[which(str_detect(res_df$id, ",", negate = FALSE)),]$id)
  # Deals with multiple ids again but at higher level
  n = 1
  number_ids = max(str_count(res_df$id, pattern = ","),na.rm = T) + 1
  if(!is.na(any(str_detect(res_df$id, ",", negate = FALSE)))){
    while (n<=number_ids & any(str_detect(res_df$id, ",", negate = FALSE))) {
      df_temp = df[which(str_detect(res_df$id, ",", negate = FALSE)),]
      res_df_temp = res_df[which(str_detect(res_df$id,",", negate = FALSE)),]

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
  }
  res_df[ranks] = getTaxonomy(res_df$id, sqlFile, desiredTaxa = str_to_lower(ranks))
  res_df = res_df[,colnames(res_df) %in% c(ranks,non_taxo_ranks)]
  res_df$superkingdom = res_df$superkingdom %>% replace_na("Unknown")
  res_df = res_df %>% column_to_rownames("feature_id")
  res_df=cbind(res_df,df[,colnames(df) %in% non_taxo_ranks,drop = FALSE])

  if(keepSAR & ("kingdom" %in% ranks)){
    index = which(Reduce(`|`, lapply(df[-1], grepl, pattern="SAR")))
    res_df[index,]$kingdom = "SAR"
  }

  if("phyloseq" %in% class(obj)){
    obj@tax_table@.Data = res_df %>% as.matrix()
    return(obj)
  }else{
    return(res_df %>% rownames_to_column("feature_id"))
  }
}
