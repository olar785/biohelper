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


blastn_taxo_assignment <- function(
    blastapp_path,
    method = "megablast",
    queries,
    megablast_opts = "-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8",
    blastn_opts = "-evalue 0.001 -max_target_seqs 5 -perc_identity 0.7",
    db,
    output_path,
    nthreads,
    minSim = 97,
    minCov = 80,
    update = FALSE,
    pident = "before",
    pgenus = 95,
    pfamily = 87,
    porder = 83,
    pclass = 81,
    pphylum = 79,
    pkingdom = 71,
    taxonly = "TRUE"
) {
  if (!dir.exists(output_path)) dir.create(output_path)

  args_base <- paste(
    "-db", db,
    "-query", queries,
    "-num_threads", nthreads
  )

  blast_outfmt <- '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'

  # Run BLAST searches
  if (method == "both" || method == "megablast") {
    message("Running megablast...")
    megablast_out <- file.path(output_path, "megablast_output.csv")
    system2(blastapp_path, args = c(args_base, megablast_opts, "-task megablast", "-out", megablast_out, blast_outfmt))
  }
  if (method == "both" || method == "blastn") {
    message("Running blastn...")
    blastn_out <- file.path(output_path, "blastn_output.csv")
    system2(blastapp_path, args = c(args_base, blastn_opts, "-task blastn", "-out", blastn_out, blast_outfmt))
  }

  # Process results using R-native LCA function
  if (method == "both") {
    lcaPident(
      blast_file = megablast_out,
      output_file = file.path(output_path, "megablast_output_processed.csv"),
      minSim = minSim, minCov = minCov, pident = pident,
      pgenus = pgenus, pfamily = pfamily, porder = porder,
      pclass = pclass, pphylum = pphylum, pkingdom = pkingdom,
      taxonly = taxonly, update = update, verbose = TRUE
    )
    lcaPident(
      blast_file = blastn_out,
      output_file = file.path(output_path, "blastn_output_processed.csv"),
      minSim = minSim, minCov = minCov, pident = pident,
      pgenus = pgenus, pfamily = pfamily, porder = porder,
      pclass = pclass, pphylum = pphylum, pkingdom = pkingdom,
      taxonly = taxonly, update = update, verbose = TRUE
    )
    message("Taxonomic assignment for both methods completed.")
  } else if (method == "megablast") {
    lcaPident(
      blast_file = megablast_out,
      output_file = file.path(output_path, "megablast_output_processed.csv"),
      minSim = minSim, minCov = minCov, pident = pident,
      pgenus = pgenus, pfamily = pfamily, porder = porder,
      pclass = pclass, pphylum = pphylum, pkingdom = pkingdom,
      taxonly = taxonly, update = update, verbose = TRUE
    )
    message("Taxonomic assignment with megablast completed.")
  } else if (method == "blastn") {
    lcaPident(
      blast_file = blastn_out,
      output_file = file.path(output_path, "blastn_output_processed.csv"),
      minSim = minSim, minCov = minCov, pident = pident,
      pgenus = pgenus, pfamily = pfamily, porder = porder,
      pclass = pclass, pphylum = pphylum, pkingdom = pkingdom,
      taxonly = taxonly, update = update, verbose = TRUE
    )
    message("Taxonomic assignment with blastn completed.")
  } else {
    stop("Invalid method. Choose from 'megablast', 'blastn', or 'both'.")
  }

  # Merging results from blastn and megablast
  if(method =="both"){
    cat("\nMerging results from blast output\n")
    megablast = fread(paste0(output_path,"/megablast_output_processed.csv")) %>% as.data.frame() %>% dplyr::mutate(method = "megablast", colsplit(taxonomy,";", names = ranks))
    blastn = fread(paste0(output_path,"/blastn_output_processed.csv")) %>% as.data.frame() %>% dplyr::mutate(method = "blastn",colsplit(taxonomy,";", names = ranks))

    blast = rbind(megablast, blastn) %>% dplyr::select(- c(taxonomy))
    blast$nRb = rowSums(is.na(blast[,ranks]) | blast[,ranks] == "")

    blast[blast==""]=NA

    newdf = data.frame(matrix(nrow=blast$ASVs %>% unique() %>% length(), ncol = length(c("ASV", ranks))))
    colnames(newdf) = c("ASV", ranks)
    newdf$ASV = blast$ASVs %>% unique()

    pb = txtProgressBar(min = 0, max = nrow(newdf), initial = 0,  style = 3)
    for (i in 1:nrow(newdf)) {
      temp_blast = blast %>% dplyr::filter(ASVs == newdf$ASV[i])
      if(temp_blast$method %>% unique() %>% length() >1){
        x = temp_blast %>% summarise(across(where(is.character), \(x) n_distinct(x, na.rm = T))) %>% dplyr::select(-c(ASVs,method))
        if(all(x %in% c(0,1))){
          index_tax = 4+length(ranks) # temp_blast initially contains 4 columns (ASVs, taxonomy, Percent_Identity, Sequence_coverage) + taxonomy split by ';'
          newdf[i,2:length(colnames(newdf))] = temp_blast[which.min(temp_blast$nRb),5:index_tax]
        }else{
          max_shared_rank = length(ranks) - max(temp_blast$nRb)
          newdf[i,2:c(max_shared_rank+1)] = temp_blast[which.min(temp_blast$nRb),5:(max_shared_rank+4)]
        }
      }
      setTxtProgressBar(pb,i)
      close(pb)
    }

  }else if(method == "megablast"){
    newdf = fread(paste0(output_path,"/megablast_output_processed.csv")) %>%
      as.data.frame() %>%
      dplyr::mutate(colsplit(taxonomy,";", names = ranks)) %>%
      dplyr::select(c("ASVs", all_of(ranks))) %>%
      dplyr::mutate(dplyr::across(where(~ is.logical(.x) & all(is.na(.x))), as.character)) %>%
      dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(.x, "")))
  }else{
    newdf = fread(paste0(output_path,"/blastn_output_processed.csv")) %>%
      as.data.frame() %>%
      dplyr::mutate(colsplit(taxonomy,";", names = ranks)) %>%
      dplyr::select(c("ASVs", all_of(ranks))) %>%
      dplyr::mutate(dplyr::across(where(~ is.logical(.x) & all(is.na(.x))), as.character)) %>%
      dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(.x, "")))
  }

  newdf$nR = rowSums(!is.na(newdf[,ranks]))
  temp_summary = newdf %>% dplyr::summarise(mean = round(mean(nR), 2), sd = round(sd(nR), 2))
  cat("\nMean assigned taxonomic ranks: ", temp_summary$mean %>%
        as.numeric(), "\nStandard deviation: ", temp_summary$sd %>%
        as.numeric(), "\n\n")
  write.table(x = newdf %>% dplyr::select(-nR), file = paste0(output_path, "/blastn_taxo_assingment.csv"), row.names = F)
  #return(newdf %>% dplyr::select(-nR))
}
