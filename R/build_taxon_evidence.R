#' Build a standard taxon evidence table
#'
#' @description
#' Creates a standardised taxon evidence table that can be passed to
#' `flag_taxa(taxon_evidence = ...)`. This scaffold does not query WoRMS, OBIS,
#' GBIF, BOLD, OpenAI, or any external service. It only extracts taxonomy,
#' chooses the most informative query taxon for each row, creates optional
#' placeholder evidence rows, validates user-supplied evidence, and returns the
#' combined evidence table.
#'
#' @param x A taxonomy table as a data frame, matrix, phyloseq taxonomy table,
#'   or a phyloseq/speedyseq object containing a taxonomy table.
#' @param tax_ranks Character vector of taxonomic ranks used to choose the query
#'   taxon, ordered from most specific to broadest.
#' @param evidence_sources Character vector naming evidence sources requested
#'   for future lookup. These sources are stored in placeholder rows but are not
#'   queried by this scaffold.
#' @param user_evidence Optional user-supplied evidence table in the same
#'   standard format accepted by `flag_taxa(taxon_evidence = ...)`.
#' @param include_empty Logical scalar. If `TRUE`, return one placeholder row
#'   for each unique query taxon even when no evidence has been fetched yet.
#' @param verbose Logical scalar. If `TRUE`, emit a short message about the
#'   number of placeholder rows created.
#'
#' @return A data frame in the standard taxon evidence format.
#'
#' @section Taxon evidence table:
#' Required columns are: `taxon_name`, `taxon_rank`, `source`,
#' `evidence_type`, `evidence_summary`, and `reference`.
#'
#' Standard optional columns included in the output are: `accepted_name`,
#' `accepted_rank`, `source_taxon_id`, `source_record_id`, `environment`,
#' `habitat`, `region`, `locality`, `decimal_latitude`, `decimal_longitude`,
#' `basis_of_record`, `occurrence_count`, `reference_url`, `doi`, and
#' `checked_at`.
#'
#' @export
#'
#' @examples
#' tax <- data.frame(
#'   kingdom = "Animalia",
#'   phylum = "Chordata",
#'   family = "Salmonidae",
#'   genus = "Salmo",
#'   species = "Salmo salar"
#' )
#'
#' evidence <- build_taxon_evidence(tax)
build_taxon_evidence <- function(
  x,
  tax_ranks = c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom"
  ),
  evidence_sources = c("worms", "obis", "gbif", "bold", "literature"),
  user_evidence = NULL,
  include_empty = TRUE,
  verbose = FALSE
) {
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")
  evidence_sources <- .validate_character_vector(
    evidence_sources,
    "evidence_sources"
  )
  user_evidence <- validate_taxon_evidence(user_evidence)
  include_empty <- .validate_logical_scalar(include_empty, "include_empty")
  verbose <- .validate_logical_scalar(verbose, "verbose")

  tax_table <- extract_tax_table(x)
  placeholder_evidence <- .empty_taxon_evidence_table()

  if (isTRUE(include_empty)) {
    placeholder_evidence <- .build_placeholder_taxon_evidence(
      tax_table = tax_table,
      tax_ranks = tax_ranks,
      evidence_sources = evidence_sources
    )
  }

  if (isTRUE(verbose)) {
    message(
      "Created ",
      nrow(placeholder_evidence),
      " placeholder taxon evidence row(s)."
    )
  }

  evidence <- rbind(
    .standardise_taxon_evidence_columns(placeholder_evidence),
    .standardise_taxon_evidence_columns(user_evidence)
  )

  validate_taxon_evidence(evidence)
}

.build_placeholder_taxon_evidence <- function(tax_table, tax_ranks, evidence_sources) {
  query_taxa <- lapply(seq_len(nrow(tax_table)), function(row_index) {
    choose_query_taxon(tax_table[row_index, , drop = FALSE], tax_ranks = tax_ranks)
  })

  query_taxa <- data.frame(
    taxon_name = vapply(
      query_taxa,
      function(query_taxon) query_taxon$query_name,
      character(1)
    ),
    taxon_rank = vapply(
      query_taxa,
      function(query_taxon) query_taxon$query_rank,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  query_taxa <- query_taxa[
    vapply(query_taxa$taxon_name, .is_informative_taxon, logical(1)) &
      vapply(query_taxa$taxon_rank, .is_informative_taxon, logical(1)),
    ,
    drop = FALSE
  ]
  query_taxa <- unique(query_taxa)

  if (nrow(query_taxa) == 0) {
    return(.empty_taxon_evidence_table())
  }

  placeholder <- .empty_taxon_evidence_table(nrow(query_taxa))
  placeholder$taxon_name <- query_taxa$taxon_name
  placeholder$taxon_rank <- query_taxa$taxon_rank
  placeholder$source <- paste(evidence_sources, collapse = ";")
  placeholder$evidence_type <- "not_queried"
  placeholder$evidence_summary <- "No evidence fetched yet."
  placeholder$reference <- NA_character_
  placeholder$checked_at <- as.character(Sys.Date())

  placeholder
}

.empty_taxon_evidence_table <- function(n = 0) {
  columns <- .taxon_evidence_standard_columns()
  out <- as.data.frame(
    stats::setNames(
      replicate(length(columns), rep(NA_character_, n), simplify = FALSE),
      columns
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out
}

.standardise_taxon_evidence_columns <- function(taxon_evidence) {
  if (is.null(taxon_evidence)) {
    return(.empty_taxon_evidence_table())
  }

  out <- .empty_taxon_evidence_table(nrow(taxon_evidence))
  shared_columns <- intersect(colnames(taxon_evidence), colnames(out))
  out[shared_columns] <- lapply(taxon_evidence[shared_columns], as.character)
  out
}

.taxon_evidence_standard_columns <- function() {
  c(
    .taxon_evidence_required_columns(),
    .taxon_evidence_optional_columns()
  )
}
