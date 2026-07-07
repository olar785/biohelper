#' Extract taxa for evidence retrieval
#'
#' @description
#' Extracts taxon names from a taxonomy table or phyloseq/speedyseq object so
#' they can be queried against evidence sources such as WoRMS. This helper only
#' prepares query taxa. It does not call WoRMS, standardise taxon names, or make
#' retain/review/exclude decisions.
#'
#' By default, the function returns unique taxa for evidence querying: the
#' highest available taxonomic assignment is chosen for each feature and
#' duplicate `taxon_name`/`taxon_rank` pairs are removed. Use
#' `include_feature_id = TRUE` and `unique = FALSE` to keep a feature-level
#' mapping table with one highest-resolution assignment per ASV/feature. Use
#' `mode = "all"` to return every valid assignment across the selected ranks.
#'
#' `extract_taxa_for_evidence()` is intentionally separate from
#' `fetch_worms_evidence()`: first extract the taxa you want to query, then pass
#' the resulting table to an evidence retrieval function.
#'
#' @param x A taxonomy table as a data frame or matrix, a phyloseq taxonomy
#'   table, or a phyloseq/speedyseq object containing a taxonomy table.
#' @param ranks Character vector of taxonomy ranks to extract. Rank columns are
#'   matched case-insensitively. The default excludes `Domain`,
#'   `Superkingdom`, and `Kingdom` because those ranks are usually too broad for
#'   evidence querying. They can still be used if explicitly included.
#' @param mode Extraction mode. `"highest"` returns the most specific valid
#'   assignment per feature. `"all"` returns all valid assignments across
#'   `ranks`.
#' @param include_feature_id Logical scalar. If `TRUE`, include a `feature_id`
#'   column using an existing feature ID column where available, otherwise using
#'   non-default taxonomy row names.
#' @param drop_uncultured Logical scalar. If `TRUE`, remove uninformative values
#'   such as `"uncultured"`, `"unidentified"`, `"unclassified"`, `"unknown"`,
#'   `"environmental sample"`, `"metagenome"`, `"Incertae sedis"`, and empty
#'   rank prefixes such as `"g__"` or `"s__"`.
#' @param drop_empty Logical scalar. If `TRUE`, remove missing, blank, and
#'   literal `"NA"` values.
#' @param unique Logical scalar. If `TRUE`, deduplicate rows by `taxon_name` and
#'   `taxon_rank`. When `feature_id` is present, the first feature ID for each
#'   duplicate taxon/rank pair is retained; `feature_id` is not part of the
#'   deduplication key.
#'
#' @return A tibble with `taxon_name` and `taxon_rank`, and optionally
#'   `feature_id`.
#' @export
#'
#' @examples
#' tax_table <- data.frame(
#'   Kingdom = "Animalia",
#'   Phylum = "Chordata",
#'   Class = "Actinopteri",
#'   Family = "Salmonidae",
#'   Genus = "Salmo",
#'   Species = "Salmo salar"
#' )
#'
#' taxa_to_query <- extract_taxa_for_evidence(tax_table)
#'
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
extract_taxa_for_evidence <- function(
  x,
  ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  mode = c("highest", "all"),
  include_feature_id = FALSE,
  drop_uncultured = TRUE,
  drop_empty = TRUE,
  unique = TRUE
) {
  ranks <- .validate_character_vector(ranks, "ranks")
  mode <- match.arg(mode)
  include_feature_id <- .validate_logical_scalar(include_feature_id, "include_feature_id")
  drop_uncultured <- .validate_logical_scalar(drop_uncultured, "drop_uncultured")
  drop_empty <- .validate_logical_scalar(drop_empty, "drop_empty")
  unique <- .validate_logical_scalar(unique, "unique")

  tax_table <- extract_tax_table(x)
  rank_matches <- .match_evidence_rank_columns(tax_table, ranks)
  if (nrow(rank_matches) == 0) {
    return(.empty_taxa_for_evidence(include_feature_id = include_feature_id))
  }

  feature_id <- .taxa_for_evidence_feature_id(tax_table)
  if (identical(mode, "highest")) {
    return(.extract_highest_taxa_for_evidence(
      tax_table = tax_table,
      rank_matches = rank_matches,
      feature_id = feature_id,
      include_feature_id = include_feature_id,
      drop_empty = drop_empty,
      drop_uncultured = drop_uncultured,
      unique = unique
    ))
  }

  .extract_all_taxa_for_evidence(
    tax_table = tax_table,
    rank_matches = rank_matches,
    feature_id = feature_id,
    include_feature_id = include_feature_id,
    drop_empty = drop_empty,
    drop_uncultured = drop_uncultured,
    unique = unique
  )
}

.extract_highest_taxa_for_evidence <- function(
  tax_table,
  rank_matches,
  feature_id,
  include_feature_id,
  drop_empty,
  drop_uncultured,
  unique
) {
  rank_matches <- .order_evidence_ranks_most_specific_first(rank_matches)
  rows <- vector("list", nrow(tax_table))
  row_count <- 0

  for (feature_index in seq_len(nrow(tax_table))) {
    for (rank_index in seq_len(nrow(rank_matches))) {
      column_name <- rank_matches$column_name[[rank_index]]
      taxon_name <- as.character(tax_table[[column_name]][[feature_index]])
      if (!.is_valid_evidence_taxon(taxon_name, drop_empty, drop_uncultured)) {
        next
      }

      row_count <- row_count + 1
      rows[[row_count]] <- data.frame(
        taxon_name = taxon_name,
        taxon_rank = rank_matches$taxon_rank[[rank_index]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      if (isTRUE(include_feature_id)) {
        rows[[row_count]]$feature_id <- feature_id[[feature_index]]
      }
      break
    }
  }

  rows <- rows[seq_len(row_count)]
  if (length(rows) == 0) {
    return(.empty_taxa_for_evidence(include_feature_id = include_feature_id))
  }

  out <- do.call(rbind, rows)
  .finalise_taxa_for_evidence(
    out = out,
    include_feature_id = include_feature_id,
    unique = unique
  )
}

.extract_all_taxa_for_evidence <- function(
  tax_table,
  rank_matches,
  feature_id,
  include_feature_id,
  drop_empty,
  drop_uncultured,
  unique
) {
  rows <- lapply(seq_len(nrow(rank_matches)), function(rank_index) {
    column_name <- rank_matches$column_name[[rank_index]]
    out <- data.frame(
      taxon_name = as.character(tax_table[[column_name]]),
      taxon_rank = rank_matches$taxon_rank[[rank_index]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (isTRUE(include_feature_id)) {
      out$feature_id <- feature_id
    }
    out
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  keep <- vapply(
    out$taxon_name,
    .is_valid_evidence_taxon,
    logical(1),
    drop_empty = drop_empty,
    drop_uncultured = drop_uncultured
  )
  out <- out[keep, , drop = FALSE]

  .finalise_taxa_for_evidence(
    out = out,
    include_feature_id = include_feature_id,
    unique = unique
  )
}

.finalise_taxa_for_evidence <- function(out, include_feature_id, unique) {
  if (nrow(out) == 0) {
    return(.empty_taxa_for_evidence(include_feature_id = include_feature_id))
  }
  if (isTRUE(unique)) {
    unique_key <- paste(out$taxon_name, out$taxon_rank, sep = "\r")
    out <- out[!duplicated(unique_key), , drop = FALSE]
  }
  if (isTRUE(include_feature_id)) {
    out <- out[, c("feature_id", "taxon_name", "taxon_rank"), drop = FALSE]
  } else {
    out <- out[, c("taxon_name", "taxon_rank"), drop = FALSE]
  }
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

.order_evidence_ranks_most_specific_first <- function(rank_matches) {
  rank_order <- .evidence_rank_specificity_order()
  specificity <- match(tolower(rank_matches$taxon_rank), tolower(rank_order))
  rank_matches[order(specificity, decreasing = TRUE), , drop = FALSE]
}

.match_evidence_rank_columns <- function(tax_table, ranks) {
  canonical_ranks <- .evidence_canonical_ranks()
  requested <- .canonicalise_evidence_ranks(ranks)
  requested <- requested[!is.na(requested)]
  if (length(requested) == 0) {
    return(data.frame(column_name = character(), taxon_rank = character()))
  }

  column_lookup <- stats::setNames(colnames(tax_table), tolower(colnames(tax_table)))
  rows <- lapply(requested, function(rank) {
    lookup_name <- tolower(rank)
    if (!(lookup_name %in% names(column_lookup))) {
      return(NULL)
    }
    column_name <- column_lookup[[lookup_name]]
    data.frame(
      column_name = column_name,
      taxon_rank = canonical_ranks[[tolower(rank)]],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(column_name = character(), taxon_rank = character()))
  }

  out <- do.call(rbind, rows)
  out[!duplicated(tolower(out$column_name)), , drop = FALSE]
}

.canonicalise_evidence_ranks <- function(ranks) {
  canonical <- .evidence_canonical_ranks()
  unname(canonical[tolower(ranks)])
}

.evidence_canonical_ranks <- function() {
  stats::setNames(
    .evidence_rank_specificity_order(),
    tolower(.evidence_rank_specificity_order())
  )
}

.evidence_rank_specificity_order <- function() {
  c(
    "Domain",
    "Superkingdom",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
}

.taxa_for_evidence_feature_id <- function(tax_table) {
  feature_columns <- c(
    "feature_id",
    "feature.id",
    "featureid",
    "taxon_id",
    "asv",
    "otu"
  )
  matches <- match(feature_columns, tolower(colnames(tax_table)))
  matches <- matches[!is.na(matches)]
  if (length(matches) > 0) {
    return(as.character(tax_table[[matches[[1]]]]))
  }

  rep(NA_character_, nrow(tax_table))
}

.is_empty_evidence_taxon <- function(taxon_name) {
  value <- trimws(as.character(taxon_name))
  is.na(taxon_name) | !nzchar(value) | tolower(value) %in% c("na", "n/a", "nan", "null")
}

.is_valid_evidence_taxon <- function(taxon_name, drop_empty, drop_uncultured) {
  if (isTRUE(drop_empty) && .is_empty_evidence_taxon(taxon_name)) {
    return(FALSE)
  }
  if (isTRUE(drop_uncultured) && .is_uncultured_evidence_taxon(taxon_name)) {
    return(FALSE)
  }
  TRUE
}

.is_uncultured_evidence_taxon <- function(taxon_name) {
  missing_value <- is.na(taxon_name)
  value <- trimws(as.character(taxon_name))
  value[missing_value] <- ""
  lower <- tolower(value)
  empty_rank_prefix <- grepl("^[a-z]__\\s*$", lower)
  asv_or_otu_placeholder <- grepl("^(asv|otu)[_-]?[0-9]+$", lower)
  exact_placeholder <- lower %in% c(
    "uncultured",
    "unidentified",
    "unclassified",
    "unknown",
    "na",
    "n/a",
    "nan",
    "null",
    "asv",
    "otu",
    "environmental sample",
    "metagenome",
    "incertae sedis"
  )
  descriptive_placeholder <- grepl(
    "\\b(uncultured|unidentified|unclassified|unknown|environmental sample|metagenome|incertae sedis)\\b",
    lower
  )

  missing_value | empty_rank_prefix | asv_or_otu_placeholder | exact_placeholder | descriptive_placeholder
}

.empty_taxa_for_evidence <- function(include_feature_id = FALSE) {
  if (isTRUE(include_feature_id)) {
    return(tibble::as_tibble(data.frame(
      feature_id = character(),
      taxon_name = character(),
      taxon_rank = character(),
      stringsAsFactors = FALSE
    )))
  }

  tibble::as_tibble(data.frame(
    taxon_name = character(),
    taxon_rank = character(),
    stringsAsFactors = FALSE
  ))
}
