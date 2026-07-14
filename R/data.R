#' Example phyloseq object
#'
#' A small `phyloseq` object used in examples and tests.
#'
#' @format A `phyloseq` object.
#' @source Included with `biohelper` for package examples.
"ps_test_data"

#' Small eukaryotic phyloseq test object
#'
#' A small `phyloseq` object containing eukaryotic amplicon data and taxonomy.
#' This dataset is intended for package examples and tests that need a compact
#' phyloseq object with eukaryotic taxonomic ranks.
#'
#' @format A `phyloseq` object with 567 taxa and 33 samples.
#'
#' @details
#' `ps_test_data_euk` is useful for testing taxonomy extraction, phyloseq-style
#' workflows, and on-demand taxon evidence preparation without relying on a
#' large external dataset.
#'
#' @examples
#' data("ps_test_data_euk")
#' taxa_to_query <- extract_taxa_for_evidence(ps_test_data_euk)
#'
#' \dontrun{
#' worms_evidence <- fetch_worms_evidence(
#'   taxa_to_query,
#'   by = "name",
#'   cache_path = "worms_cache.rds"
#' )
#' }
#'
#' @source Internal package test data.
#' @keywords datasets
"ps_test_data_euk"
