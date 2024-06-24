#' Loads taxo_normalisation
#'
#' @description
#' This function performs normalisation of taxonomic assignments of eukaryotic data
#' identified via silva, midori or any other reference databases, using NCBI curated
#' classification and nomenclature database. It can also provide further taxonomic
#' information of additional desired ranks. It can either take in a taxonomy table
#' (e.g. from Qiime2) or a phyloseq object. For species assignation, make sure that
#' the species column contains both the Genus and Species names. If it doesn't, use "spcn = T".
#'
#' Importantly, this function requires the download of NCBI taxonomic database.
#'
#' This takes few minutes to install using the following command:\cr
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package. The database is ~65 GB but the user can set getAccessions=FALSE to drastically reduce its size.
#'
#' @param
#' obj = Either a dataframe containing taxonomic information (e.g. output of Qiime2 or the assignTaxonomy from dada2) or a phyloseq object
#' @param
#' sqlFile = Path to the local NCBI taxonomy db
#' @param
#' ranks = Ranks to return
#' @param
#' keepSAR = Keep the SAR assignment from the input data, which is not a valid group in NCBI taxonomy db. If false, SAR taxa are kept but the 'SAR' assignation is removed. Note that some databases may not use the 'SAR' label so if merging data assigned with different databases and KeepSAR is TRUE, there may be discrepancies in the data (default is FALSE)
#' @param
#' spnc = Only needs to be applied if the Genus is not present under the Species column (default is FALSE).
#' @export
#' @examples
#' data("ps_test_data")
#' taxo_normalisation(ps_test_data, sqlFile = 'accessionTaxa.sql')
#'

taxo_normalisation = function(obj, sqlFile, keepSAR = F, spnc = F, ranks = c("Superkingdom", "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus", "Species")){

  if("phyloseq" %in% class(obj)){
    df = obj@tax_table@.Data %>% as.data.frame()
  } else{
    df = obj %>% as.data.frame()
  }

  df <- df[,!is.na(colnames(df))]
  ranks = tolower(ranks)
  colnames(df) = tolower(colnames(df))
  taxa_id = c("otu","otus","asv","asvs","feature_id","feature.id","nR")

  # Formatting the df
  if("taxon" %in% colnames(df)){
    df$taxon = gsub("D_\\d+__","",taxo$taxon)
    df = df %>% dplyr::select(c(colnames(df)[which(colnames(df)%in%c(taxa_id,"taxon"))])) %>% splitstackshape::cSplit('taxon', sep = ';')
  }else if("taxonomy" %in% colnames(df)){
    df = df %>% dplyr::select(c(colnames(df)[which(colnames(df)%in%c(taxa_id,"taxonomy"))])) %>% splitstackshape::cSplit('taxonomy', sep = ';')
  }
  df = df %>% dplyr::mutate_if(is.factor, as.character)
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

  paternsToRemove = c("^.+_environmental.+|environmental_.+|uncultured_.+|_sp\\..+|_sp.|_sp.+| sp\\..+| sp.| sp.+|_\\(.+|^.+_metagenome|_cf.|s__|g__|f__|c__|o__|c__|p__|'")
  df = df %>% dplyr::mutate_all(list(~str_replace(.,paternsToRemove, ""))) %>% dplyr::mutate_all(list(~na_if(.,"")))
  if("species" %in% colnames(df)){
    if(spnc == F){
      df = df %>% dplyr::mutate(species = paste0(sapply(strsplit(species,"_"), `[`, 1)," ", sapply(strsplit(species,"_"), `[`, 2)),
                                species = gsub("NA","",species),
                                species = gsub("^\\s+$",NA,species)) %>%
        dplyr::mutate_all(list(~na_if(.,"")))
    }
    if(spnc == T){
      df = df %>% dplyr::mutate(species = paste0(genus," ",species),
                                species = gsub("NA","",species),
                                species = gsub("^\\s+$",NA,species),
                                species = gsub("^\\s+(.+)","\\1",species)) %>%
        dplyr::mutate_all(list(~na_if(.,"")))
    }
  }

  df = df %>% dplyr::mutate_all(list(~str_replace(.,"NA NA| NA", ""))) %>% dplyr::mutate(across(everything(), gsub, pattern = "_", replacement = " ")) %>% dplyr::mutate_all(list(~na_if(.,""))) %>%
    dplyr::mutate(across(everything(), gsub, pattern = "^\\s+|\\s+$", replacement = "")) # The latter is to removing leading and trailing spaces

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
    res_df[which(is.na(res_df$id)),]$id = ifelse(!is.na(id),id,NA)
    r = r + 1
  }
  # Deals with multiple ids
  n = 1
  number_ids = max(str_count(res_df$id, pattern = ","),na.rm = T) + 1
  while (n<=number_ids & any(str_detect(res_df$id, ",", negate = FALSE), na.rm = T)) {
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
  if(!is.na(any(str_detect(res_df$id, ",", negate = FALSE), na.rm = T))){
    while (n<=number_ids & any(str_detect(res_df$id, ",", negate = FALSE), na.rm = T)) {
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
    if(length(index)>0){
      res_df[index,]$kingdom = "SAR"
    }
  }

  if("phyloseq" %in% class(obj)){
    obj@tax_table@.Data = res_df %>% as.matrix()
    return(obj)
  }else{
    return(res_df %>% rownames_to_column("feature_id"))
  }
}
