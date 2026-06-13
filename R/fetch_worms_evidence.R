#' Fetch WoRMS taxon evidence for selected taxa
#'
#' @description
#' Queries WoRMS for a supplied list of taxon names or AphiaIDs and returns the
#' standard taxon evidence format used by `flag_taxa()` and
#' `build_taxon_evidence()`. This function does not download or build the full
#' WoRMS database. It is intended for enriching a small set of taxa and for
#' building a reusable local WoRMS evidence cache.
#'
#' For larger enrichment runs, query gently, use `sleep` where repeated AphiaID
#' requests are needed, and provide `cache_path` so successful results can be
#' reused across sessions.
#'
#' @param taxa Character vector of taxon names when `by = "name"`, or an
#'   integer, numeric, or character vector of AphiaIDs when `by = "aphia_id"`.
#' @param by Query mode. Use `"name"` to query taxon names with
#'   `worrms::wm_records_names()`, or `"aphia_id"` to query AphiaIDs with
#'   `worrms::wm_record()`.
#' @param marine_only Logical scalar passed to `worrms::wm_records_names()` when
#'   `by = "name"`.
#' @param cache_path Optional path to an RDS cache. Existing cached rows are
#'   reused and new query results are appended.
#' @param sleep Numeric scalar. Seconds to pause between repeated AphiaID
#'   requests.
#' @param checked_at Date or character scalar recorded in the `checked_at`
#'   column.
#' @param verbose Logical scalar. If `TRUE`, emit short progress messages.
#'
#' @return A data frame in the standard taxon evidence format.
#'
#' @details
#' The `worrms` package is a suggested dependency. If it is not installed and a
#' live query is needed, `fetch_worms_evidence()` errors with an installation
#' message. Tests and examples should mock WoRMS responses or use cached data
#' rather than making live API calls.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' evidence <- fetch_worms_evidence(
#'   taxa = c("Salmo salar", "Gadus morhua"),
#'   by = "name",
#'   cache_path = "worms_evidence_cache.rds"
#' )
#' }
fetch_worms_evidence <- function(
  taxa,
  by = c("name", "aphia_id"),
  marine_only = FALSE,
  cache_path = NULL,
  sleep = 0.2,
  checked_at = Sys.Date(),
  verbose = TRUE
) {
  by <- match.arg(by)
  taxa <- .validate_worms_taxa(taxa, by)
  marine_only <- .validate_logical_scalar(marine_only, "marine_only")
  cache_path <- .validate_optional_cache_path(cache_path)
  sleep <- .validate_non_negative_numeric_scalar(sleep, "sleep")
  checked_at <- .validate_checked_at(checked_at)
  verbose <- .validate_logical_scalar(verbose, "verbose")

  if (length(taxa) == 0) {
    return(.empty_taxon_evidence_table())
  }

  cache <- .read_worms_evidence_cache(cache_path)
  cached_indices <- .match_worms_cache(taxa, by, cache)
  taxa_to_query <- unique(taxa[is.na(cached_indices)])

  new_evidence <- .empty_taxon_evidence_table()
  if (length(taxa_to_query) > 0) {
    .require_worrms()
    if (isTRUE(verbose)) {
      message("Querying WoRMS for ", length(taxa_to_query), " taxon/taxa.")
    }
    new_evidence <- .query_worms_evidence(
      taxa = taxa_to_query,
      by = by,
      marine_only = marine_only,
      sleep = sleep,
      checked_at = checked_at
    )
    cache <- rbind(cache, new_evidence)
    .write_worms_evidence_cache(cache_path, cache)
    cached_indices <- .match_worms_cache(taxa, by, cache)
  }

  out <- cache[cached_indices, , drop = FALSE]
  rownames(out) <- NULL
  validate_taxon_evidence(out)
}

.require_worrms <- function() {
  if (!requireNamespace("worrms", quietly = TRUE)) {
    stop(
      "Live WoRMS queries require the `worrms` package. ",
      "Install it with `install.packages(\"worrms\")`, or provide a ",
      "`cache_path` containing the requested taxa.",
      call. = FALSE
    )
  }
}

.worms_records_names <- function(taxa, marine_only = FALSE) {
  worrms::wm_records_names(taxa, marine_only = marine_only)
}

.worms_record <- function(aphia_id) {
  if (!exists("wm_record", envir = asNamespace("worrms"), inherits = FALSE)) {
    stop(
      "The installed `worrms` package does not provide `wm_record()` for ",
      "AphiaID lookups.",
      call. = FALSE
    )
  }

  worrms::wm_record(aphia_id)
}

.query_worms_evidence <- function(
  taxa,
  by,
  marine_only,
  sleep,
  checked_at
) {
  if (identical(by, "name")) {
    return(.query_worms_evidence_by_name(
      taxa = taxa,
      marine_only = marine_only,
      checked_at = checked_at
    ))
  }

  .query_worms_evidence_by_aphia_id(
    taxa = taxa,
    sleep = sleep,
    checked_at = checked_at
  )
}

.query_worms_evidence_by_name <- function(taxa, marine_only, checked_at) {
  records <- tryCatch(
    .worms_records_names(taxa, marine_only = marine_only),
    error = function(error) {
      stop(
        "WoRMS name query failed: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  rows <- lapply(seq_along(taxa), function(query_index) {
    query <- taxa[[query_index]]
    candidates <- .worms_records_for_name_query(
      records = records,
      query = query,
      query_index = query_index,
      query_count = length(taxa)
    )
    selected <- .select_worms_name_record(candidates, query)
    if (is.null(selected$record)) {
      return(.worms_no_match_evidence(query, by = "name", checked_at = checked_at))
    }

    .worms_record_to_evidence(
      record = selected$record,
      query = query,
      by = "name",
      match_type = selected$match_type,
      checked_at = checked_at
    )
  })

  validate_taxon_evidence(do.call(rbind, rows))
}

.query_worms_evidence_by_aphia_id <- function(taxa, sleep, checked_at) {
  rows <- lapply(seq_along(taxa), function(query_index) {
    if (query_index > 1 && sleep > 0) {
      Sys.sleep(sleep)
    }

    query <- taxa[[query_index]]
    record <- tryCatch(
      .worms_record(query),
      error = function(error) {
        stop(
          "WoRMS AphiaID query failed for `",
          query,
          "`: ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    record <- .as_worms_records_data_frame(record)
    if (nrow(record) == 0) {
      return(.worms_no_match_evidence(query, by = "aphia_id", checked_at = checked_at))
    }

    .worms_record_to_evidence(
      record = record[1, , drop = FALSE],
      query = query,
      by = "aphia_id",
      match_type = "aphia_id",
      checked_at = checked_at
    )
  })

  validate_taxon_evidence(do.call(rbind, rows))
}

.worms_records_for_name_query <- function(records, query, query_index, query_count) {
  if (is.null(records)) {
    return(data.frame())
  }

  if (!inherits(records, "data.frame") && is.list(records)) {
    record_names <- names(records)
    if (.is_worms_record_like_list(records)) {
      return(.as_worms_records_data_frame(records))
    }
    if (!is.null(record_names)) {
      named_index <- which(tolower(record_names) == tolower(query))
      if (length(named_index) > 0) {
        return(.as_worms_records_data_frame(records[[named_index[[1]]]]))
      }
    }
    if (length(records) == query_count || length(records) == 1) {
      return(.as_worms_records_data_frame(records[[query_index]]))
    }

    return(.as_worms_records_data_frame(records))
  }

  records <- .as_worms_records_data_frame(records)
  if (nrow(records) == 0) {
    return(records)
  }
  if (query_count == 1) {
    return(records)
  }

  query_columns <- c(
    "query",
    "query_name",
    "input",
    "input_name",
    "submitted_name",
    "name_submitted",
    "requested_name"
  )
  query_values <- .row_first_non_empty(records, query_columns)
  has_query_values <- .is_non_empty_value(query_values)
  if (any(has_query_values)) {
    matched_query <- tolower(trimws(query_values)) == tolower(trimws(query))
    return(records[matched_query, , drop = FALSE])
  }

  exact <- .worms_exact_name_match(records, query)
  records[exact, , drop = FALSE]
}

.select_worms_name_record <- function(records, query) {
  records <- .as_worms_records_data_frame(records)
  if (nrow(records) == 0) {
    return(list(record = NULL, match_type = "no_match"))
  }

  exact_accepted <- .worms_exact_name_match(records, query) &
    .worms_accepted_status(records)
  if (any(exact_accepted)) {
    return(list(
      record = records[which(exact_accepted)[[1]], , drop = FALSE],
      match_type = "accepted_exact"
    ))
  }

  list(
    record = records[1, , drop = FALSE],
    match_type = if (isTRUE(.worms_exact_name_match(records[1, , drop = FALSE], query))) {
      "exact"
    } else {
      "first_returned"
    }
  )
}

.worms_record_to_evidence <- function(record, query, by, match_type, checked_at) {
  record <- .as_worms_records_data_frame(record)
  environment <- .worms_environment_summary(record)
  source_taxon_id <- .worms_record_value(record, c("AphiaID", "aphia_id"))
  source_record_id <- .worms_record_value(record, c("valid_AphiaID", "valid_aphia_id"))
  matched_name <- .worms_record_value(record, c(
    "scientificname",
    "scientificName",
    "valid_name",
    "accepted_name"
  ))
  accepted_name <- .worms_record_value(record, c(
    "valid_name",
    "accepted_name",
    "scientificname",
    "scientificName"
  ))
  rank <- .worms_record_value(record, c("rank", "taxonRank"))
  status <- .worms_record_value(record, c("status", "valid_status"))
  reference <- .worms_record_value(record, c("citation", "reference"))

  if (!.is_non_empty_value(matched_name) && identical(by, "aphia_id")) {
    matched_name <- query
  }
  if (!.is_non_empty_value(rank)) {
    rank <- "unknown"
  }
  if (!.is_non_empty_value(status)) {
    status <- "unknown"
  }
  if (!.is_non_empty_value(reference)) {
    reference <- "WoRMS"
  }

  evidence <- .empty_taxon_evidence_table(1)
  evidence$taxon_name <- if (identical(by, "name")) query else matched_name
  evidence$taxon_rank <- rank
  evidence$source <- "worms"
  evidence$evidence_type <- "taxonomic_environment_database"
  evidence$evidence_summary <- .worms_evidence_summary(
    query = query,
    matched_name = matched_name,
    status = status,
    match_type = match_type,
    environment = environment
  )
  evidence$reference <- reference
  evidence$accepted_name <- accepted_name
  evidence$accepted_rank <- rank
  evidence$source_taxon_id <- source_taxon_id
  evidence$source_record_id <- source_record_id
  evidence$environment <- environment
  evidence$reference_url <- .worms_reference_url(record, source_taxon_id)
  evidence$checked_at <- checked_at

  evidence
}

.worms_no_match_evidence <- function(query, by, checked_at) {
  evidence <- .empty_taxon_evidence_table(1)
  evidence$taxon_name <- query
  evidence$taxon_rank <- "unknown"
  evidence$source <- "worms"
  evidence$evidence_type <- "taxonomic_environment_database"
  evidence$evidence_summary <- "No WoRMS record found for queried taxon."
  evidence$reference <- "WoRMS"
  if (identical(by, "aphia_id")) {
    evidence$source_taxon_id <- query
  }
  evidence$checked_at <- checked_at

  evidence
}

.worms_evidence_summary <- function(
  query,
  matched_name,
  status,
  match_type,
  environment
) {
  environment_text <- if (.is_non_empty_value(environment)) {
    paste0("; environment flags: ", environment)
  } else {
    ""
  }

  paste0(
    "Queried taxon '",
    query,
    "'; matched WoRMS taxon '",
    matched_name,
    "'; status: ",
    status,
    "; match_type: ",
    match_type,
    environment_text,
    "."
  )
}

.as_worms_records_data_frame <- function(records) {
  if (is.null(records)) {
    return(data.frame())
  }
  if (inherits(records, "data.frame")) {
    return(as.data.frame(records, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (is.list(records) && length(records) == 0) {
    return(data.frame())
  }
  if (is.list(records) && all(vapply(records, inherits, logical(1), "data.frame"))) {
    return(do.call(rbind, lapply(records, as.data.frame, stringsAsFactors = FALSE)))
  }

  as.data.frame(records, stringsAsFactors = FALSE, check.names = FALSE)
}

.is_worms_record_like_list <- function(records) {
  record_names <- names(records)
  if (is.null(record_names)) {
    return(FALSE)
  }

  any(tolower(record_names) %in% tolower(c(
    "AphiaID",
    "scientificname",
    "scientificName",
    "valid_name",
    "rank",
    "status"
  )))
}

.worms_exact_name_match <- function(records, query) {
  candidate_names <- cbind(
    .source_column(records, "scientificname"),
    .source_column(records, "scientificName"),
    .source_column(records, "valid_name"),
    .source_column(records, "accepted_name")
  )

  apply(candidate_names, 1, function(names) {
    any(tolower(trimws(names)) == tolower(trimws(query)), na.rm = TRUE)
  })
}

.worms_accepted_status <- function(records) {
  status <- tolower(.row_first_non_empty(records, c("status", "valid_status")))
  status %in% c("accepted", "valid")
}

.worms_record_value <- function(record, candidates) {
  value <- .row_first_non_empty(record, candidates)
  if (length(value) == 0 || !.is_non_empty_value(value[[1]])) {
    return(NA_character_)
  }

  value[[1]]
}

.read_worms_evidence_cache <- function(cache_path) {
  if (is.null(cache_path) || !file.exists(cache_path)) {
    return(.empty_taxon_evidence_table())
  }

  cache <- readRDS(cache_path)
  if (!inherits(cache, "data.frame")) {
    stop("`cache_path` must contain a data.frame saved as RDS.", call. = FALSE)
  }
  cache <- validate_taxon_evidence(cache)
  .standardise_taxon_evidence_columns(cache)
}

.write_worms_evidence_cache <- function(cache_path, cache) {
  if (is.null(cache_path)) {
    return(invisible(cache))
  }

  cache <- .standardise_taxon_evidence_columns(cache)
  cache <- validate_taxon_evidence(cache)
  cache_dir <- dirname(cache_path)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(cache, cache_path)
  invisible(cache)
}

.match_worms_cache <- function(taxa, by, cache) {
  if (nrow(cache) == 0) {
    return(rep(NA_integer_, length(taxa)))
  }

  cache_source <- as.character(cache$source) == "worms"
  cache_type <- as.character(cache$evidence_type) == "taxonomic_environment_database"
  cache_available <- cache_source & cache_type

  vapply(taxa, function(query) {
    if (identical(by, "name")) {
      matches <- which(
        cache_available &
          tolower(trimws(as.character(cache$taxon_name))) == tolower(trimws(query))
      )
    } else {
      matches <- which(
        cache_available &
          trimws(as.character(cache$source_taxon_id)) == trimws(query)
      )
    }

    if (length(matches) == 0) {
      return(NA_integer_)
    }
    matches[[1]]
  }, integer(1))
}

.validate_worms_taxa <- function(taxa, by) {
  if (missing(taxa) || is.null(taxa)) {
    stop("`taxa` must be supplied.", call. = FALSE)
  }
  if (length(taxa) == 0) {
    return(character())
  }
  if (identical(by, "name") && !is.character(taxa)) {
    stop("`taxa` must be a character vector when `by = \"name\"`.", call. = FALSE)
  }
  if (
    identical(by, "aphia_id") &&
      !(is.character(taxa) || is.numeric(taxa) || is.integer(taxa))
  ) {
    stop(
      "`taxa` must be an integer, numeric, or character vector when ",
      "`by = \"aphia_id\"`.",
      call. = FALSE
    )
  }

  taxa <- trimws(as.character(taxa))
  if (any(is.na(taxa)) || any(!nzchar(taxa))) {
    stop("`taxa` must not contain missing or empty values.", call. = FALSE)
  }

  taxa
}

.validate_optional_cache_path <- function(cache_path) {
  if (is.null(cache_path)) {
    return(NULL)
  }
  if (
    !is.character(cache_path) ||
      length(cache_path) != 1 ||
      is.na(cache_path) ||
      !nzchar(trimws(cache_path))
  ) {
    stop("`cache_path` must be NULL or a non-empty character scalar.", call. = FALSE)
  }

  cache_path
}

.validate_non_negative_numeric_scalar <- function(value, name) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      value < 0
  ) {
    stop("`", name, "` must be a non-negative numeric scalar.", call. = FALSE)
  }

  value
}

.validate_checked_at <- function(checked_at) {
  if (
    length(checked_at) != 1 ||
      is.na(checked_at) ||
      !nzchar(trimws(as.character(checked_at)))
  ) {
    stop("`checked_at` must be a non-empty scalar.", call. = FALSE)
  }

  as.character(checked_at)
}
