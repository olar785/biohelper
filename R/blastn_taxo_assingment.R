#' blastn_taxo_assignment
#'
#' @description
#' This wrapper function takes in a fasta file and returns a dataframe of
#' taxonomic assignments. It first perform a blastn and/or megablast
#' search in a database such as NCBI. It then uses a Last Common Ancestor (LCA)
#' approach and optionally, a minimum percent identity (pident; per taxonomic rank)
#' to assign taxonomy.
#'
#' If blastn and megablast were both chosen, it merges the results
#' using the highest taxonomic resolution between the two, that is, if there is consensus across
#' all assigned ranks. Otherwise, it assigns taxonomy using the LCA approach.
#'
#' The function will write files from each step and return a dataframe of
#' taxonomic assignments.
#'
#' Megablast and especially blastn can take several hours
#' to run so if executed on a remote computer, best to execute as a sbatch job
#' or using a terminal multiplexer like screen or tmux.
#'
#' Users only interested in performing a LCA/pident taxonomic
#' assignment from multiple blast hits may refer to the lcaPident function (can be
#' particularly useful to test overall taxonomic resolution using different parameters).
#'
#' @param
#' blastapp_path      Path to the blastn program
#' @param
#' method             Method to be used for blast search. Options are 'megablast', 'blastn' and 'both' (default).
#' @param
#' queries            Fasta file to query
#' @param
#' megablast_opts     megablast options including evalue (default = 0.001), max_target_seqs (default = 5) and perc_identity (default = 0.8)
#' @param
#' blastn_opts        blastn options including evalue (default = 0.001), max_target_seqs (default = 5) and perc_identity (default = 0.5)
#' @param
#' db                 Reference database
#' @param
#' output_path        Path to output directory
#' @param
#' nthreads           Number of threads for the blast search
#' @param
#' minSim             Minimum similarity to assign species hits (default: 97)
#' @param
#' minCov             Minimum coverage to keep hits (default: 80)
#' @param
#' update             Should the taxonomy database be updated? (default: FALSE)
#' @param
#' pident             To reduce taxonomy assignment according to default percent identity thresholds. Options are: before or after LCA assingment
#' @param
#' pgenus             Minimum similarity to assign genus (default: 95)
#' @param
#' pfamily            Minimum similarity to assign family (default: 87)
#' @param
#' porder             Minimum similarity to assign order (default: 83)
#' @param
#' pclass             Minimum similarity to assign class (default: 81)
#' @param
#' pphylum            Minimum similarity to assign phylum (default: 79)
#' @param
#' pkingdom           Minimum similarity to assign kingdom (default: 71)
#'
#' @export
#' @examples
#' blastn_taxo_assignment(blastapp_path = "/home/miniconda3/bin/blastn",\cr
#' queries="uniqueSeqsTest.fasta",\cr
#' megablast_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8",\cr
#' blastn_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.5",db="nt",\cr
#' output_path="blast_results",\cr
#' nthreads=10)


blastn_taxo_assignment = function(blastapp_path,
                            method="both",
                            queries,
                            megablast_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8",
                            blastn_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.7",
                            db,
                            output_path,
                            nthreads,
                            minSim=97,
                            minCov=80,
                            update=FALSE,
                            pident="no",
                            pgenus=95,
                            pfamily=87,
                            porder=83,
                            pclass=81,
                            pphylum=79,
                            pkingdom=71,
                            taxonly="TRUE")
  {
  # Blastn and Megablast
  if(!dir.exists(output_path)) dir.create(output_path)
  args <- paste(paste("-db", db, collapse = " "),
                paste("-query", queries, collapse = " "),
                paste("-num_threads", nthreads, collapse = " "))

  if(method == "both"){
    cat("\nPerforming megablast and blastn\n")
    system2(blastapp_path, args = c(args,megablast_opts,"-task megablast", paste0("-out"," ",output_path,"/megablast_output.csv"), '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'))
    cat("\nDone with megablast, now starting blastn\n")
    system2(blastapp_path, args = c(args,blastn_opts,"-task blastn", paste0("-out"," ",output_path,"/blastn_output.csv"), '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'))
    cat("\nDone with blastn\n")
  }else if(method == "megablast"){
    cat("\nPerforming megablast only\n")
    system2(blastapp_path, args = c(args,megablast_opts,"-task megablast", paste0("-out"," ",output_path,"/megablast_output.csv"), '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'))

  }else if(method == "blastn"){
    cat("\nPerforming blastn only\n")
    system2(blastapp_path, args = c(args,blastn_opts,"-task blastn", paste0("-out"," ",output_path,"/blastn_output.csv"), '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'))

  }else{
    stop("The method specified is not available. Options are 'megablast', 'blastn' and 'both'", call. = FALSE)
  }

  # R packages used
  packages <- c("tidyverse", "data.table","reshape2","reticulate")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))

  # Taxonomic assignment using LCA and pident
  for (i in c("pandas","ete3","csv","re","argparse","numpy", "time", "tqdm")) {
    if(reticulate::py_module_available(!i)){
      reticulate::py_install(i)
    }
  }

  reticulate::py_run_file(file = system.file("env.py",package = "biohelper"))
  pyscript = system.file("Pident_LCA_blast_taxo_assignment.py",package = "biohelper")

  args_general = paste(
  paste("--minSim", minSim, collapse = " "),
  paste("--minCov", minCov, collapse = " "),
  paste("--update", update, collapse = " "),
  paste("--pident", pident, collapse = " "),
  paste("--pgenus", pgenus, collapse = " "),
  paste("--pfamily", pfamily, collapse = " "),
  paste("--porder", porder, collapse = " "),
  paste("--pclass", pclass, collapse = " "),
  paste("--pphylum", pphylum, collapse = " "),
  paste("--pkingdom", pkingdom, collapse = " "),
  paste("--taxonly", taxonly, collapse = " "))

  cat("\nPerforming taxonomic assignment from blast output\n")
  args_megablast = paste(
    paste("-b", paste0(output_path,"/megablast_output.csv"), collapse = " "),
    paste("-o", paste0(output_path,"/megablast_output_processed.csv"), collapse = " "),
    args_general)
  args_blastn = paste(
    paste("-b", paste0(output_path,"/blastn_output.csv"), collapse = " "),
    paste("-o", paste0(output_path,"/blastn_output_processed.csv"), collapse = " "),
    args_general)

  if(method == "both"){
    system2(pyscript, args_megablast)
    system2(pyscript, args_blastn)
  }else if(method == "megablast"){
    system2(pyscript, args_megablast)
    newdf = fread(paste0(output_path,"/megablast_output_processed.csv")) %>% dplyr::mutate(colsplit(taxonomy,";", names = c("superkingdom","kingdom","phylum","class","order","family","genus","species"))) %>% dplyr::select(-c(taxonomy,Percent_Identity,Sequence_coverage)) %>% mutate_all(na_if,"")
  }else{
    system2(pyscript, args_blastn)
    newdf = fread(paste0(output_path,"/blastn_output_processed.csv")) %>% dplyr::mutate(colsplit(taxonomy,";", names = c("superkingdom","kingdom","phylum","class","order","family","genus","species"))) %>% dplyr::select(-c(taxonomy,Percent_Identity,Sequence_coverage)) %>% mutate_all(na_if,"")
  }

  # Merging results from blastn and megablast
  if(method =="both"){
    megablast = fread(paste0(output_path,"/megablast_output_processed.csv")) %>% dplyr::mutate(method = "megablast", colsplit(taxonomy,";", names = c("superkingdom","kingdom","phylum","class","order","family","genus","species")))
    blastn = fread(paste0(output_path,"/blastn_output_processed.csv")) %>% dplyr::mutate(method = "blastn",colsplit(taxonomy,";", names = c("superkingdom","kingdom","phylum","class","order","family","genus","species")))

    blast = rbind(megablast, blastn) %>% dplyr::select(- c(taxonomy))
    blast$nRb = rowSums(is.na(blast[,c("superkingdom","kingdom","phylum","class","order","family","genus","species")] ) | blast[,c("superkingdom","kingdom","phylum","class","order","family","genus","species")] == "")
    blast <- blast %>% mutate_all(na_if,"")

    newdf = data.frame(matrix(nrow=blast$ASVs %>% unique() %>% length(), ncol = 9))
    colnames(newdf) = c("ASV","superkingdom","kingdom","phylum","class","order","family","genus","species")
    newdf$ASV = blast$ASVs %>% unique()

    for (i in 1:nrow(newdf)) {
      temp_blast = blast %>% dplyr::filter(ASVs == newdf$ASV[i])
      if(temp_blast$method %>% unique() %>% length() >1){
        x = temp_blast %>% summarise(across(where(is.character), n_distinct,na.rm = T)) %>% dplyr::select(-c(ASVs,method))
        if(all(x %in% c(0,1))){
          newdf[i,2:9] = temp_blast[which.min(temp_blast$nRb),5:12]
        }else{
          max_shared_rank = 8 - max(temp_blast$nRb)
          newdf[i,2:c(max_shared_rank+1)] = temp_blast[which.min(temp_blast$nRb),5:(max_shared_rank+4)]
        }
      }
    }
  }
  newdf$nRb = rowSums(is.na(newdf[,c("superkingdom","kingdom","phylum","class","order","family","genus","species")] ) | newdf[,c("superkingdom","kingdom","phylum","class","order","family","genus","species")] == "")
  temp_summary = newdf %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric())
  write.table(x = newdf %>% dplyr::select(-nRb), file = paste0(output_path,"/blastn_taxo_assingment.csv"))
  return(newdf %>% dplyr::select(-nRb))
}



