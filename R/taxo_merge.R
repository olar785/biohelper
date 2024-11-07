#' taxo_merge
#'
#' @description
#' This function takes in a list of taxonomic assignment dataframes, normalises the
#' taxonomy using NCBI curated taxonomic database, and uses the highest
#' taxonomic resolution between the assignments if there is consensus across
#' all assigned ranks and dataframes. Otherwise, it assigns taxonomy the Last common
#' ancestor (LCA) if there is a consensus (>50%) across dataframes.
#'
#' Importantly, this function requires the download of NCBI taxonomic database (~65 GB).
#' Additionally, the first column of the dataframe needs to be the feature id (e.g. ASV, OTU, MOTU, ZOTU identifier, etc.).
#'
#' This takes few minutes to install using the following command:
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package.
#'
#' @param
#' df_list        List of dataframes with taxonomic assignment
#' @param
#' sqlFile        Path to NCBI taxonomic reference database
#' @param
#' ranks          Ranks to return
#' @param
#' addExtra       Currently adds the Protozoa and Archaeplastida (excl. Viridiplantae) groups under the Kingdom rank (default is TRUE). This is to make it easier to dissociate the Plantea and Protozoa group from Metazoans
#'
#' @export
#' @examples
#' taxo_merge(df_list = list(df1, df2),  sqlFile = path_to_NCBI_taxo_db, ranks = c("Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"), addExtra=T)


taxo_merge = function(
    df_list,
    sqlFile,
    ranks,
    addExtra = T)
{
  ranks = tolower(ranks)

  for (i in 1:length(df_list)) {
    cat("\nNormalizing taxonomy of object #",i,"\n")
    colnames(df_list[[i]]) = colnames(df_list[[i]]) %>% tolower()
    colnames(df_list[[i]])[1] = "feature_id"
    df_list[[i]] = taxo_normalisation(obj = df_list[[i]], sqlFile = sqlFile, ranks = ranks, addExtra = addExtra)
    df_list[[i]] = df_list[[i]] %>%
      dplyr::mutate(df=as.character(i)) %>%
      dplyr::mutate_all(list(~na_if(., "Unknown")))
  }

  dfall = data.table::rbindlist(df_list)
  dfall$nRb = rowSums(is.na(dfall[,..ranks] ) | dfall[,..ranks] == "")
  #dfall <- dfall %>% mutate_all(na_if,"")
  dfall <- dfall %>% dplyr::mutate(across(where(is.character), ~ na_if(.x, "")))
  dfall = dfall[!nRb==length(ranks),]

  newdf = data.frame(matrix(nrow=dfall$feature_id %>% unique() %>% length(), ncol = length(ranks)+1))
  colnames(newdf) = c("feature_id",ranks)
  newdf$feature_id = dfall$feature_id %>% unique()

  pb = txtProgressBar(min = 0, max = nrow(newdf), initial = 0,  style = 3)
  for (i in 1:nrow(newdf)) {
    temp = dfall %>% dplyr::filter(feature_id == newdf$feature_id[i])
    if(temp$df %>% unique() %>% length() >1){
      x = temp %>% dplyr::summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
      while (!all(x %in% c(0,1))) {
        r = min(which(apply(temp, 2, function(x) length(unique(x[!is.na(x)]))) > 1))
        if(which.max(table(temp[,..r])) > nrow(temp)/2){
          temp = temp[which(temp[,..r] == names(which.max(table(temp[,..r])))),]
          x = temp %>% dplyr::summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
        }else{
          temp[,r:which(colnames(temp)==tail(ranks,1))] = NA
          x = temp %>% dplyr::summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
        }
      }
      newdf[i,1:(length(ranks)+1)] = temp[which.min(temp$nRb),1:(length(ranks)+1)]
    }else{
      newdf[i,1:(length(ranks)+1)] = temp[which.min(temp$nRb),1:(length(ranks)+1)]
    }
    setTxtProgressBar(pb,i)
    close(pb)
  }

  newdf$nRb =  rowSums(!is.na(newdf[,ranks] ) & newdf[,ranks] != "")
  temp_summary = newdf %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric(),"\n\n")
  return(newdf %>% dplyr::select(-nRb))
}
