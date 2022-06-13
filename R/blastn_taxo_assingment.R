#' Loads blastn_taxo_assignment
#'
#' This function takes in a fasta file and returns taxonomic assignments.
#' Specifically, it performs blastn and megablast in the nt database of NCBI and
#' then uses a Last Common Ancestor (LCA) approach and a minimum percent identity
#' (per taxonomic rank) to assign taxonomy. It then merges the results from blastn
#' and megablast using the highest taxonomic resolution with a consensus between
#' the two approaches.
#'
#' @export
#' @examples
#' blast_assignment(blastapp_path = "/home/miniconda3/bin/blastn",
#' queries="uniqueSeqs.fasta",
#' megablast_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8",
#' blastn_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.5",db="nt",
#' output="blast_resutls",
#' nthreads=20)


blastn_taxo_assignment = function(blastapp_path,
                            queries,
                            megablast_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8",
                            blastn_opts="-evalue 0.001 -max_target_seqs 5 -perc_identity 0.5",
                            db,
                            output_path,
                            nthreads,
                            minSim=97,
                            minCov=80,
                            update=FALSE,
                            pident="no",
                            taxonly="TRUE")
  {
  # Blastn and Megablast
  if(!dir.exists(output_path)) dir.create(output_path)
  args <- paste(paste("-db", db, collapse = " "),
                paste("-query", queries, collapse = " "),
                paste("-num_threads", nthreads, collapse = " "))
  system2(blastapp_path, args = c(args,megablast_opts,"-task megablast", paste0("-out"," ",output_path,"/megablast_output.csv"), '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'))
  system2(blastapp_path, args = c(args,blastn_opts,"-task blastn", paste0("-out"," ",output_path,"/blastn_output.csv"), '-outfmt "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"'))

  # Taxonomic assignment using LCA and pident
  pyscript = system.file("Pident_LCA_blast_taxo_assignment.py",package = "biohelper")

  args_general = paste(
  paste("--minSim", minSim, collapse = " "),
  paste("--minCov", minCov, collapse = " "),
  paste("--update", update, collapse = " "),
  paste("--pident", pident, collapse = " "),
  paste("--taxonly", taxonly, collapse = " "))

  args_megablast = paste(
    paste("-b", paste0(output_path,"/megablast_output.csv"), collapse = " "),
    paste("-o", paste0(output_path,"/megablast_output_processed.csv"), collapse = " "),
    args_general
  )
  args_blastn = paste(
    paste("-b", paste0(output_path,"/blastn_output.csv"), collapse = " "),
    paste("-o", paste0(output_path,"/blastn_output_processed.csv"), collapse = " "),
    args_general
  )

  system2(pyscript, args_megablast)
  system2(pyscript, args_blastn)

  # Merging results from blastn and megablast
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
  return(newdf)
}











