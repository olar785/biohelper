#' taxo_merge
#'
#' @description
#' This function takes in two taxonomic assignment dataframes, normalises the
#' taxonomy using NCBI curated taxonomic database, and uses the highest
#' taxonomic resolution between the two assignments if there is consensus across
#' all assigned ranks. Otherwise, it assigns taxonomy using a Last common
#' ancestor approach (LCA).
#'
#' Importantly, this function requires the download of NCBI taxonomic database (~65 GB).
#'
#' This takes few minutes to install using the following command:
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package.
#'
#' @param
#' df1            First dataframe with taxonomic assignment
#' @param
#' df2            Second dataframe with taxonomic assignment
#' @param
#' sqlFile        Path to NCBI taxonomic reference database
#' @param
#' ranks          Ranks to return
#' @param
#' keepSAR        Keep the SAR assignment from the input data, which is not a valid group in NCBI taxonomy db (default is FALSE)
#'
#' @export
#' @examples
#' taxo_merge(df1 = df1, df2 = df2,  sqlFile = path_to_NCBI_taxo_db, ranks = c("Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"), keepSAR=F)


taxo_merge = function(df1,
                      df2,
                      sqlFile,
                      ranks,
                      keepSAR = F)
{
  ranks = tolower(ranks)

  colnames(df1) = colnames(df1) %>% tolower()
  colnames(df1)[1] = "feature_id"
  df1 = taxo_normalisation(obj = df1, sqlFile = sqlFile, ranks = ranks,keepSAR = keepSAR)

  colnames(df2) = colnames(df2) %>% tolower()
  colnames(df2)[1] = "feature_id"
  df2 = taxo_normalisation(obj = df2, sqlFile = sqlFile, ranks = ranks, keepSAR = keepSAR)

  dfall = rbind(df1 %>% dplyr::mutate(df='1'), df2 %>% dplyr::mutate(df='2'))
  dfall$nRb = rowSums(is.na(dfall[,ranks] ) | dfall[,ranks] == "")
  dfall <- dfall %>% mutate_all(na_if,"")

  newdf = data.frame(matrix(nrow=dfall$feature_id %>% unique() %>% length(), ncol = length(ranks)+1))
  colnames(newdf) = c("feature_id",ranks)
  newdf$feature_id = dfall$feature_id %>% unique()

  for (i in 1:nrow(newdf)) {
    temp = dfall %>% dplyr::filter(feature_id == newdf$feature_id[i])
    if(temp$df %>% unique() %>% length() >1){
      x = temp %>% summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
      if(all(x %in% c(0,1))){
        newdf[i,1:length(ranks)+1] = temp[which.min(temp$nRb),1:length(ranks)+1]
      }else{
        max_shared_rank = length(ranks) - max(temp$nRb)
        newdf[i,1:c(max_shared_rank+1)] = temp[which.min(temp$nRb),1:(max_shared_rank+1)]
      }
    }
  }

  newdf$nRb = rowSums(!is.na(newdf[,ranks] ) & newdf[,ranks] != "")
  temp_summary = newdf %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric())
  return(newdf %>% dplyr::select(-nRb))
}
