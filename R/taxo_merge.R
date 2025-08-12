#' taxo_merge
#'
#' @description
#' This function takes in a list of taxonomic assignment dataframes, normalizes the
#' taxonomy using NCBI curated taxonomic database, and uses the highest
#' taxonomic resolution between the assignments if there is consensus across
#' all assigned ranks and dataframes. Otherwise, it assigns taxonomy the Last common
#' ancestor (LCA) if there is a consensus (>50%) across dataframes.
#' By using the 'priority_df' option, it is also possible to use the taxonomy of a specific dataframe (assignment) instead of using LCA.
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
#' priority_df    Taxonomic dataframe to prioritize in case of discrepancy between assignments (default is NA).
#' @param
#' addExtra       Currently adds the Protozoa and Archaeplastida (excl. Viridiplantae) groups under the Kingdom rank (default is TRUE). This is to make it easier to dissociate the Plantea and Protozoa group from Metazoans
#' @param
#' spnc           Only needs to be applied if the Genus is not present under the Species column (default is FALSE).
#'
#' @export
#' @examples
#' taxo_merge(df_list = list(df1, df2),  sqlFile = path_to_NCBI_taxo_db, ranks = c("Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"), priority_df = NA, addExtra=T, spnc = F)


taxo_merge = function(
    df_list,
    sqlFile,
    ranks = c("Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    priority_df = NA,
    addExtra = T,
    spnc = F)
{
  ranks = tolower(ranks)

  for (i in 1:length(df_list)) {
    cat("\nNormalizing taxonomy of object #",i,"\n")
    colnames(df_list[[i]]) = colnames(df_list[[i]]) %>% tolower()
    colnames(df_list[[i]])[1] = "feature_id"
    df_list[[i]] = taxo_normalisation(obj = df_list[[i]], sqlFile = sqlFile, ranks = ranks, addExtra = addExtra, spnc = F)
    df_list[[i]] = df_list[[i]] %>%
      dplyr::mutate(df=as.character(i)) %>%
      dplyr::mutate_all(list(~na_if(., "Unknown")))
  }

  dfall = data.table::rbindlist(df_list) %>% as.data.frame()
  dfall$nRb = rowSums(is.na(dfall[,ranks] ) | dfall[,ranks] == "")
  #dfall <- dfall %>% mutate_all(na_if,"")
  dfall <- dfall %>% dplyr::mutate(across(where(is.character), ~ na_if(.x, "")))
  dfall = dfall[!dfall$nRb==length(ranks),]

  newdf = data.frame(matrix(nrow=dfall$feature_id %>% unique() %>% length(), ncol = length(ranks)+1))
  colnames(newdf) = c("feature_id",ranks)
  newdf$feature_id = dfall$feature_id %>% unique()

  cat("\nMerging taxonomy\n")
  pb = txtProgressBar(min = 0, max = nrow(newdf), initial = 0,  style = 3)
  for (i in 1:nrow(newdf)) {
    # df per ASV
    temp = dfall %>% dplyr::filter(feature_id == newdf$feature_id[i]) %>% as.data.table()
    # If more than 1 taxo assignment result for this ASV...
    if(temp$df %>% unique() %>% length() >1){
      # If a specific df should be used in case of discrepancy...
      if(!is.na(priority_df)){
        # Assign taxonomy from that df
        temp = temp[temp$df==as.character(priority_df)]
        newdf[i,1:(length(ranks)+1)] = temp[which.min(temp$nRb),1:(length(ranks)+1)]
      }else{
        # Determine the number of different values at each rank
        x = temp %>% dplyr::summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
        # While some ranks have more than 1 different values...
        while (!all(x %in% c(0,1))) {
          # Determine the lowest resolution rank where there is discrepancy
          r = min(which(apply(temp, 2, function(x) length(unique(x[!is.na(x)]))) > 1))
          # If one of the assignments at that rank has majority...
          if(which.max(table(temp[,..r])) > nrow(temp)/2){
            # Assign the complete taxonomy from the assignment that has majority
            temp = temp[which(temp[,..r] == names(which.max(table(temp[,..r])))),]
            # Determine the number of different values at each rank
            x = temp %>% dplyr::summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
          }else{
            # If no assignment has majority, replace that assignment by NA
            temp[,r:which(colnames(temp)==tail(ranks,1))] = NA
            # Determine the number of different values at each rank
            x = temp %>% dplyr::summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(feature_id,df))
          }
        }
      }
      # Once consensus is reached, assigned taxonomy from the assignment with highest resolution
      newdf[i,1:(length(ranks)+1)] = temp[which.min(temp$nRb),1:(length(ranks)+1)]
      # In case there is already consensus across all assignment, keep the one with highest resolution
    }else{
      newdf[i,1:(length(ranks)+1)] = temp[which.min(temp$nRb),1:(length(ranks)+1)]
    }

    setTxtProgressBar(pb,i)
    close(pb)
  }
  # Determine the number of rank assigned for each ASV
  newdf$nRb =  rowSums(!is.na(newdf[,ranks] ) & newdf[,ranks] != "")
  # Provide the mean and sd of number of rank assigned across the consensus taxo data
  temp_summary = newdf %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric(),"\n\n")
  return(newdf %>% dplyr::select(-nRb))
}
