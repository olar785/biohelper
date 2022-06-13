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
                            update=TRUE,
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

  b_megablast = paste0(output_path,"/megablast_output.csv")
  b_blastn = paste0(output_path,"/blastn_output.csv")
  o_megablast = paste0(output_path,"/megablast_output_processed.csv")
  o_blastn = paste0(output_path,"/blastn_output_processed.csv")

  args_megablast = paste(
    paste("-b", b_megablast, collapse = " "),
    paste("-o", o_megablast, collapse = " "),
    paste("--minSim", minSim, collapse = " "),
    paste("--minCov", minCov, collapse = " "),
    paste("--update", update, collapse = " "),
    paste("--pident", pident, collapse = " "),
    paste("--taxonly", taxonly, collapse = " ")
  )
  system2(pyscript, args_megablast)

  args_blastn = paste(
    paste("-b", b_blastn, collapse = " "),
    paste("-o", o_blastn, collapse = " "),
    paste("--minSim", minSim, collapse = " "),
    paste("--minCov", minCov, collapse = " "),
    paste("--update", update, collapse = " "),
    paste("--pident", pident, collapse = " "),
    paste("--taxonly", taxonly, collapse = " ")
  )
  system2(pyscript, args_blastn)

  # Merging results from blastn and megablast
  blastn = fread("")
  megablast = fread("")


  }











