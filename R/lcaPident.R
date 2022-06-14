#' Loads lcaPident
#'
#' @description
#' This  function takes in the output of a blastn search with multiple hits and
#' uses a Last Common Ancestor (LCA) approach and optionally, a minimum percent
#' identity (pident; per taxonomic rank) to assign taxonomy.
#'
#' @param
#' blastn_file        Blast file containing multiple hits per query sequence. To work with this function the blastn search MUST have the format output 6  with options "query.id", "query.length", "pident", "subject.id", "subject.GBid", "evalue", "bit.score","staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp"
#' minSim             Minimum similarity to assign species hits (default: 97)
#' minCov             Minimum coverage to keep hits (default: 80)
#' update             Should the taxonomy database be updated? (default: FALSE)
#' pident             To reduce taxonomy assingment according to default percent identity thresholds. Options are: before or after LCA assingment
#' pgenus             Minimum similarity to assign genus (default: 95)
#' pfamily            Minimum similarity to assign family (default: 87)
#' porder             Minimum similarity to assign order (default: 83)
#' pclass             Minimum similarity to assign class (default: 81)
#' pphylum            Minimum similarity to assign phylum (default: 79)
#' pkingdom           Minimum similarity to assign kingdom (default: 71)
#'
#' @export
#' @examples
#' lcaPident(blastn_file = blastn_file, output = "taxo_assignment.csv", pident="no")


lcaPident = function(blastn_file,
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

use_condaenv(condaenv = "r-reticulate")

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
  paste("-b", blastn_file, collapse = " "),
  paste("-o", output, collapse = " "),
  args_general
)
system2(pyscript, args_blastn)
temp = fread(output) %>% dplyr::mutate(method = "blastn",colsplit(taxonomy,";", names = c("superkingdom","kingdom","phylum","class","order","family","genus","species")))
temp_summary = temp %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(sd(nRb),2))
cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric())
return(temp %>% dplyr::select(-nRb))
}
