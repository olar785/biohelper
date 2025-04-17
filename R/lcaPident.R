#' lcaPident
#'
#' @description
#' This function takes in the output of a blastn search with multiple hits and
#' uses a Last Common Ancestor (LCA) approach and optionally, a minimum percent
#' identity (pident; per taxonomic rank) to assign taxonomy.
#'
#' @param
#' blast_file        Blast file containing multiple hits per query sequence. To work with this function the blastn search MUST have the format output 6  with options "query.id", "query.length", "pident", "subject.id", "subject.GBid", "evalue", "bit.score","staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp"
#' @param
#' minSim             Minimum similarity to assign species hits (default: 97)
#' @param
#' minCov             Minimum coverage to keep hits (default: 80)
#' @param
#' update             Should the taxonomy database be updated? (default: FALSE)
#' @param
#' pident             To reduce taxonomy assingment according to default percent identity thresholds. Options are: before or after LCA assingment
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
#' lcaPident(blast_file = blastn_file, output = "taxo_assignment.csv", pident="no")


lcaPident = function(blast_file,
                     output,
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
    if(!reticulate::py_module_available(i)){
      reticulate::py_install(i)
    }
  }

  reticulate::py_run_file(file = system.file("env.py",package = "biohelper"))


  # Taxonomic assignment using LCA and pident
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

  args_blastn = paste(
    paste("-b", blast_file, collapse = " "),
    paste("-o", output, collapse = " "),
    args_general
  )
  #system2(pyscript, args_blastn)
  # Run Python script and capture exit status
  exit_status <- system2(pyscript, args_blastn, stdout = TRUE, stderr = TRUE)

  # Check for failure
  if (!is.null(attr(exit_status, "status")) && attr(exit_status, "status") != 0) {
    stop("Python script failed to execute. Aborting.")
  }

  ranks = c("domain","superkingdom","kingdom","phylum","class","order","family","genus","species")
  temp = fread(output) %>% dplyr::mutate(colsplit(taxonomy,";", names = ranks))
  temp$nRb = rowSums(!is.na(temp[, ..ranks] ) & temp[, ..ranks] != "")

  temp_summary = temp %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
  cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric())
  return(temp %>% dplyr::select(-c(taxonomy,nRb,Percent_Identity,Sequence_coverage)))
}
