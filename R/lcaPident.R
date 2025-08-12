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
#' pident             To reduce taxonomy assignment according to default percent identity thresholds. Options are: before or after LCA assignment
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
#' lcaPident(blast_file = blastn_file, output = "taxo_assignment.csv", pident="before")

# Function to assign taxonomy from BLAST file

lcaPident <- function(
    blast_file, output_file, ftbl_file = NULL,
    minSim = 97, minCov = 80, pident = "no",
    pgenus = 95, pfamily = 87, porder = 83,
    pclass = 81, pphylum = 79, pkingdom = 71,
    taxonly = TRUE, update = FALSE, verbose = FALSE
) {

  desired_ranks <- c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")

  # Optionally update the NCBI taxonomy database
  if (update) {
    message("Updating local NCBI taxonomy database...")
    taxizedb::db_download_ncbi(overwrite = TRUE)
  }
  # Load BLAST table
  # Read with fread for reliable parsing, disable quote handling to avoid warnings
  blast <- data.table::fread(
    blast_file,
    col.names = c(
      "query.id", "query.length", "pident", "subject.id", "subject.GBid",
      "evalue", "bit.score", "staxids", "sscinames", "sblastnames",
      "qcovs", "qcovhsp"
    ),
    quote = ""
  )

  # Trim leading/trailing quotes if present
  blast <- blast %>% dplyr::mutate(
    dplyr::across(everything(), ~stringr::str_remove_all(.x, '^"|"$'))
  )

  # Remove unwanted entries
  blast_trimmed <- dplyr::filter(blast, !stringr::str_detect(sscinames, regex("uncultured|unidentified|environmental sample", ignore_case = TRUE))) %>%
    dplyr::filter(qcovs >= minCov)

  # Get taxonomy for staxids using lapply
  ## Keeping the first staxid if multiple
  staxids = sub(";.*", "", blast_trimmed$staxid)
  staxid_list <- unique(staxids)
  cat("Number of unique staxid =", length(staxid_list), "\n")

  pb <- progress::progress_bar$new(total = length(staxid_list), format = "[:bar] :current/:total (:percent) eta: :eta", clear = FALSE)
  taxo_result <- lapply(staxid_list, function(tx) {
    pb$tick()
    res <- tryCatch({
      taxizedb::classification(tx, db = "ncbi")
    }, error = function(e) NULL)

    if (is.null(res) || is.null(res[[1]]) || length(res[[1]]) == 0 || all(is.na(res[[1]]))) {
      message("[!] No taxonomy found for staxid: ", tx)
      tax_vec <- stats::setNames(rep("NA", length(desired_ranks)), desired_ranks)
    } else {
      ranks_df <- dplyr::filter(res[[1]], rank %in% desired_ranks)
      tax_vec <- stats::setNames(rep("NA", length(desired_ranks)), desired_ranks)
      tax_vec[ranks_df$rank] <- ranks_df$name
    }
    dplyr::tibble(staxids = tx, !!!as.list(tax_vec))
  })
  taxo_df <- dplyr::bind_rows(taxo_result)
  blast_annotated <- dplyr::left_join(blast_trimmed, taxo_df, by = "staxids")

  # Apply pident thresholds - each threshold affects only its corresponding rank
  apply_pident_thresholds <- function(df, minSim, pk, pp, pc, po, pf, pg) {
    df <- df %>% dplyr::mutate(across(all_of(desired_ranks), as.character))
    # Each threshold affects only its corresponding rank
    df <- df %>% dplyr::mutate(
      kingdom = ifelse(pident < pk, "NA", kingdom),
      phylum = ifelse(pident < pp, "NA", phylum),
      class = ifelse(pident < pc, "NA", class),
      order = ifelse(pident < po, "NA", order),
      family = ifelse(pident < pf, "NA", family),
      genus = ifelse(pident < pg, "NA", genus),
      species = ifelse(pident < minSim, "NA", species)
    )
    return(df)
  }

  # LCA logic
  calculate_lca <- function(tbl) {
    ranks <- rev(desired_ranks)
    tbl_lca <- dplyr::group_by(tbl, query.id) %>%
      dplyr::group_map(~ {
        group_df <- .x
        group_key <- .y  # This contains the grouping variable (query.id)

        best_hit <- dplyr::slice_min(group_df, order_by = evalue, with_ties = FALSE)
        lca <- best_hit[1,]

        # Add back the query.id from the grouping key
        lca$query.id <- group_key$query.id

        for (rank in ranks) {
          taxa <- unique(group_df[[rank]])
          # Filter out NA values, empty strings, and "NA" strings
          valid_taxa <- taxa[!is.na(taxa) & taxa != "" & taxa != "NA"]

          if (length(valid_taxa) > 1) {
            # Set current rank and all downstream ranks to "NA"
            rank_index <- which(desired_ranks == rank)
            lca[desired_ranks[rank_index:length(desired_ranks)]] <- "NA"
            break
          }
        }
        lca
      }) %>%
      dplyr::bind_rows()

    tbl_lca <- tidyr::unite(tbl_lca, "taxonomy", dplyr::all_of(desired_ranks),
                            sep = ";", remove = FALSE, na.rm = TRUE) %>%
      dplyr::mutate(taxonomy = dplyr::if_else(taxonomy == "", "Unknown", taxonomy))

    return(tbl_lca)
  }

  blast_annotated$pident = as.numeric(blast_annotated$pident)

  if (pident == "before") {
    blast_annotated <- apply_pident_thresholds(blast_annotated, minSim, pk = pkingdom , pp = pphylum, pc = pclass, po = porder, pf = pfamily, pg = pgenus)
    final_tbl <- calculate_lca(blast_annotated)
  } else if (pident == "after") {
    final_tbl <- calculate_lca(blast_annotated)
    final_tbl <- apply_pident_thresholds(final_tbl, minSim, pk = pkingdom , pp = pphylum, pc = pclass, po = porder, pf = pfamily, pg = pgenus)
    final_tbl <- tidyr::unite(final_tbl, "taxonomy", dplyr::all_of(desired_ranks), sep = ";", remove = FALSE, na.rm = TRUE)
  } else {
    final_tbl <- calculate_lca(blast_annotated)
  }

  feature_table <- dplyr::select(final_tbl, query.id, taxonomy, pident, qcovs, staxids)
  colnames(feature_table) = c("ASVs", "taxonomy", "Percent_Identity", "Sequence_coverage", "staxid")
  # ASVs/OTUs with NA across all ranks either had no consensus at domain level or the staxid was not found by taxizedb

  readr::write_csv(feature_table, output_file)
  if (verbose) {
    message("\nTransformed table written to: ", output_file)
  }

  return(feature_table)
}
