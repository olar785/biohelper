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
#' Input taxa are deduplicated before cache lookup or live WoRMS querying, so
#' duplicate feature-level rows do not trigger redundant API calls. The returned
#' table contains one evidence row per unique queried taxon name or AphiaID.
#'
#' @param taxa Character vector of taxon names when `by = "name"`, a data frame
#'   containing a `taxon_name` column when `by = "name"`, or an integer,
#'   numeric, or character vector of AphiaIDs when `by = "aphia_id"`. AphiaIDs
#'   are normalised to numeric values before querying WoRMS.
#' @param by Query mode. Use `"name"` to query taxon names with
#'   `worrms::wm_records_names()`, or `"aphia_id"` to query AphiaIDs with
#'   `worrms::wm_record()`.
#' @param marine_only Logical scalar passed to `worrms::wm_records_names()` when
#'   `by = "name"`.
#' @param cache_path Optional path to an RDS cache. Existing cached rows are
#'   reused and new query results are appended.
#' @param sleep Numeric scalar. Seconds to pause between repeated AphiaID
#'   requests or live query batches.
#' @param batch_size Positive integer scalar. Maximum number of unique taxa sent
#'   to WoRMS per live query batch. Defaults to `25`.
#' @param max_tries Positive integer scalar. Maximum number of attempts for
#'   transient WoRMS/API failures per batch. Defaults to `3`.
#' @param retry_sleep Non-negative numeric scalar. Seconds to wait between retry
#'   attempts. Defaults to `5`.
#' @param continue_on_error Logical scalar. If `TRUE`, repeated batch failures
#'   return conservative failed-lookup evidence rows with
#'   `evidence_type = "taxonomy_lookup_failed"` and a warning. If `FALSE`,
#'   repeated failures abort with the original WoRMS error.
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
#' WoRMS no matches and transient lookup failures are missing evidence, not
#' evidence of ecological incompatibility.
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
  batch_size = 25,
  max_tries = 3,
  retry_sleep = 5,
  continue_on_error = TRUE,
  checked_at = Sys.Date(),
  verbose = TRUE
) {
  by <- match.arg(by)
  query_rank_lookup <- .worms_query_rank_lookup(taxa, by)
  taxa <- .validate_worms_taxa(taxa, by)
  taxa <- .deduplicate_worms_taxa(taxa, by)
  marine_only <- .validate_logical_scalar(marine_only, "marine_only")
  cache_path <- .validate_optional_cache_path(cache_path)
  sleep <- .validate_non_negative_numeric_scalar(sleep, "sleep")
  batch_size <- .validate_positive_integer_scalar(batch_size, "batch_size")
  max_tries <- .validate_positive_integer_scalar(max_tries, "max_tries")
  retry_sleep <- .validate_non_negative_numeric_scalar(retry_sleep, "retry_sleep")
  continue_on_error <- .validate_logical_scalar(continue_on_error, "continue_on_error")
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
      batch_count <- ceiling(length(taxa_to_query) / batch_size)
      message(
        "fetch_worms_evidence(): querying WoRMS for ",
        length(taxa_to_query),
        " unique taxa in ",
        batch_count,
        " batches of up to ",
        batch_size,
        "."
      )
    }
    new_evidence <- .query_worms_evidence(
      taxa = taxa_to_query,
      by = by,
      marine_only = marine_only,
      sleep = sleep,
      batch_size = batch_size,
      max_tries = max_tries,
      retry_sleep = retry_sleep,
      continue_on_error = continue_on_error,
      query_rank_lookup = query_rank_lookup,
      checked_at = checked_at,
      verbose = verbose
    )
    cacheable_evidence <- .worms_cacheable_evidence(new_evidence)
    cache <- .bind_rows_fill(list(cache, cacheable_evidence))
    .write_worms_evidence_cache(cache_path, cache)
  }

  available <- .bind_rows_fill(list(cache, new_evidence))
  evidence_indices <- .match_worms_cache(
    taxa,
    by,
    available,
    include_failed = TRUE
  )
  if (any(is.na(evidence_indices))) {
    missing <- taxa[is.na(evidence_indices)]
    stop(
      "Internal error: missing WoRMS evidence rows for queried taxa: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out <- available[evidence_indices, , drop = FALSE]
  rownames(out) <- NULL
  out <- validate_taxon_evidence(out)
  if (isTRUE(verbose)) {
    unresolved <- .worms_unresolved_evidence(out)
    message(
      "fetch_worms_evidence(): returned WoRMS evidence for ",
      nrow(out) - sum(unresolved),
      " / ",
      length(taxa),
      " unique taxa; ",
      sum(unresolved),
      " unresolved."
    )
  }
  out
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

.worms_record <- function(id, verbose = FALSE) {
  if (!exists("wm_record", envir = asNamespace("worrms"), inherits = FALSE)) {
    stop(
      "The installed `worrms` package does not provide `wm_record()` for ",
      "AphiaID lookups.",
      call. = FALSE
    )
  }

  api_ids <- .normalise_worms_aphia_ids(
    id,
    drop_invalid = FALSE,
    warn = FALSE,
    value_name = "`id`"
  )
  .call_worrms_wm_record(id = api_ids, verbose = verbose)
}

.call_worrms_wm_record <- function(id, verbose = FALSE) {
  api_ids <- id
  stopifnot(
    is.atomic(api_ids),
    !is.list(api_ids),
    is.numeric(api_ids) || is.integer(api_ids)
  )
  if (isTRUE(verbose)) {
    message(
      "Calling worrms::wm_record() with ",
      length(api_ids),
      " numeric AphiaID(s)."
    )
  }

  worrms::wm_record(id = api_ids)
}

.query_worms_evidence <- function(
  taxa,
  by,
  marine_only,
  sleep,
  batch_size,
  max_tries,
  retry_sleep,
  continue_on_error,
  query_rank_lookup,
  checked_at,
  verbose
) {
  batches <- split(taxa, ceiling(seq_along(taxa) / batch_size))
  rows <- vector("list", length(batches))
  for (batch_index in seq_along(batches)) {
    rows[[batch_index]] <- .query_worms_evidence_batch_with_retries(
      taxa = batches[[batch_index]],
      by = by,
      marine_only = marine_only,
      max_tries = max_tries,
      retry_sleep = retry_sleep,
      continue_on_error = continue_on_error,
      query_rank_lookup = query_rank_lookup,
      checked_at = checked_at,
      verbose = verbose,
      batch_index = batch_index,
      batch_count = length(batches)
    )
    if (batch_index < length(batches) && sleep > 0) {
      Sys.sleep(sleep)
    }
  }

  validate_taxon_evidence(.bind_rows_fill(rows))
}

.query_worms_evidence_batch_with_retries <- function(
  taxa,
  by,
  marine_only,
  max_tries,
  retry_sleep,
  continue_on_error,
  query_rank_lookup,
  checked_at,
  verbose,
  batch_index,
  batch_count
) {
  last_error <- NULL
  attempts_used <- 0L
  for (attempt in seq_len(max_tries)) {
    attempts_used <- attempt
    result <- tryCatch(
      .query_worms_evidence_batch(
        taxa = taxa,
        by = by,
        marine_only = marine_only,
        checked_at = checked_at,
        verbose = verbose
      ),
      error = function(error) error
    )
    if (!inherits(result, "error")) {
      return(result)
    }

    last_error <- result
    retryable <- .is_retryable_worms_error(result)
    if (isTRUE(verbose)) {
      message(
        "fetch_worms_evidence(): WoRMS batch ",
        batch_index,
        " / ",
        batch_count,
        " failed on attempt ",
        attempt,
        " / ",
        max_tries,
        ": ",
        conditionMessage(result),
        if (retryable && attempt < max_tries) " Retrying." else ""
      )
    }
    if (!retryable || attempt >= max_tries) {
      break
    }
    if (retry_sleep > 0) {
      Sys.sleep(retry_sleep)
    }
  }

  if (isTRUE(continue_on_error)) {
    warning(
      "WoRMS query batch ",
      batch_index,
      " / ",
      batch_count,
      " failed after ",
      attempts_used,
      " attempt(s): ",
      if (is.null(last_error)) "unknown error" else conditionMessage(last_error),
      ". Returning conservative failed-lookup evidence rows for this batch.",
      call. = FALSE
    )
    return(.worms_lookup_failed_evidence(
      taxa = taxa,
      by = by,
      checked_at = checked_at,
      query_rank_lookup = query_rank_lookup
    ))
  }

  stop(
    "WoRMS query batch ",
    batch_index,
    " / ",
    batch_count,
    " failed: ",
    conditionMessage(last_error),
    call. = FALSE
  )
}

.is_retryable_worms_error <- function(error) {
  text <- tolower(paste(
    conditionMessage(error),
    paste(class(error), collapse = " ")
  ))
  grepl(
    "empty reply from server|server returned nothing|timeout|timed out|connection reset|curl|network|failed to connect|could not resolve|temporar|http (500|502|503|504)|\\b(500|502|503|504)\\b",
    text
  )
}

.worms_lookup_failed_evidence <- function(
  taxa,
  by,
  checked_at,
  query_rank_lookup = character()
) {
  rows <- lapply(taxa, function(query) {
    query_key <- .worms_query_key(query, by)
    evidence <- .empty_taxon_evidence_table(1)
    evidence$taxon_name <- as.character(query)
    evidence$taxon_rank <- .worms_query_rank(query_key, query_rank_lookup)
    evidence$source <- "worms"
    evidence$evidence_type <- "taxonomy_lookup_failed"
    evidence$evidence_summary <- "WoRMS lookup failed due to transient API/network error; no compatibility inference made."
    evidence$reference <- NA_character_
    evidence$accepted_name <- NA_character_
    evidence$accepted_rank <- NA_character_
    evidence$source_taxon_id <- if (identical(by, "aphia_id")) query_key else NA_character_
    evidence$source_record_id <- NA_character_
    evidence$environment <- NA_character_
    evidence$habitat <- NA_character_
    evidence$region <- NA_character_
    evidence$locality <- NA_character_
    evidence$decimal_latitude <- NA_character_
    evidence$decimal_longitude <- NA_character_
    evidence$basis_of_record <- NA_character_
    evidence$occurrence_count <- NA_character_
    evidence$reference_url <- NA_character_
    evidence$doi <- NA_character_
    evidence$checked_at <- checked_at
    evidence
  })

  validate_taxon_evidence(.bind_rows_fill(rows))
}

.query_worms_evidence_batch <- function(
  taxa,
  by,
  marine_only,
  checked_at,
  verbose
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
    checked_at = checked_at,
    verbose = verbose
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

  validate_taxon_evidence(.bind_rows_fill(rows))
}

.query_worms_evidence_by_aphia_id <- function(taxa, checked_at, verbose) {
  api_ids <- .normalise_worms_aphia_ids(
    taxa,
    drop_invalid = TRUE,
    warn = TRUE,
    value_name = "`taxa`"
  )
  query_keys <- .worms_aphia_id_key(api_ids)

  records <- tryCatch(
    .worms_record(id = api_ids, verbose = verbose),
    error = function(error) {
      stop(
        "WoRMS AphiaID query failed for `",
        paste(query_keys, collapse = ", "),
        "`: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  records <- .as_worms_records_data_frame(records)

  rows <- lapply(seq_along(api_ids), function(query_index) {
    query_key <- query_keys[[query_index]]
    record <- .worms_records_for_aphia_id_query(
      records = records,
      query_key = query_key,
      query_index = query_index,
      query_count = length(api_ids)
    )
    if (nrow(record) == 0) {
      return(.worms_no_match_evidence(query_key, by = "aphia_id", checked_at = checked_at))
    }

    .worms_record_to_evidence(
      record = record[1, , drop = FALSE],
      query = query_key,
      by = "aphia_id",
      match_type = "aphia_id",
      checked_at = checked_at
    )
  })

  validate_taxon_evidence(.bind_rows_fill(rows))
}

.worms_records_for_aphia_id_query <- function(records, query_key, query_index, query_count) {
  records <- .as_worms_records_data_frame(records)
  if (nrow(records) == 0) {
    return(records)
  }

  record_ids <- .worms_aphia_id_key(.row_first_non_empty(records, c("AphiaID", "aphia_id")))
  valid_record_ids <- .worms_aphia_id_key(.row_first_non_empty(records, c(
    "valid_AphiaID",
    "valid_aphia_id"
  )))
  matched_id <- (!is.na(record_ids) & record_ids == query_key) |
    (!is.na(valid_record_ids) & valid_record_ids == query_key)
  if (any(matched_id)) {
    return(records[which(matched_id)[[1]], , drop = FALSE])
  }

  if (nrow(records) == query_count) {
    return(records[query_index, , drop = FALSE])
  }

  if (query_count == 1) {
    return(records[1, , drop = FALSE])
  }

  records[FALSE, , drop = FALSE]
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

.worms_cacheable_evidence <- function(evidence) {
  if (nrow(evidence) == 0) {
    return(evidence)
  }
  evidence[
    as.character(evidence$evidence_type) != "taxonomy_lookup_failed",
    ,
    drop = FALSE
  ]
}

.match_worms_cache <- function(taxa, by, cache, include_failed = FALSE) {
  if (nrow(cache) == 0) {
    return(rep(NA_integer_, length(taxa)))
  }

  cache_source <- as.character(cache$source) == "worms"
  cache_types <- "taxonomic_environment_database"
  if (isTRUE(include_failed)) {
    cache_types <- c(cache_types, "taxonomy_lookup_failed")
  }
  cache_type <- as.character(cache$evidence_type) %in% cache_types
  cache_available <- cache_source & cache_type

  if (identical(by, "aphia_id")) {
    cache_aphia_ids <- .worms_aphia_id_key(cache$source_taxon_id)
    query_aphia_ids <- .worms_aphia_id_key(taxa)
  }

  vapply(seq_along(taxa), function(query_index) {
    query <- taxa[[query_index]]
    if (identical(by, "name")) {
      matches <- which(
        cache_available &
          tolower(trimws(as.character(cache$taxon_name))) == tolower(trimws(query))
      )
    } else {
      query_key <- query_aphia_ids[[query_index]]
      matches <- which(
        cache_available &
          !is.na(cache_aphia_ids) &
          !is.na(query_key) &
          cache_aphia_ids == query_key
      )
    }

    if (length(matches) == 0) {
      return(NA_integer_)
    }
    matches[[1]]
  }, integer(1))
}

.worms_unresolved_evidence <- function(evidence) {
  if (nrow(evidence) == 0) {
    return(logical())
  }
  evidence_type <- as.character(evidence$evidence_type)
  summary <- tolower(as.character(evidence$evidence_summary))
  evidence_type == "taxonomy_lookup_failed" |
    grepl("no worms record found|lookup failed", summary)
}

.worms_query_rank_lookup <- function(taxa, by) {
  if (!inherits(taxa, "data.frame") || !identical(by, "name")) {
    return(character())
  }
  if (!all(c("taxon_name", "taxon_rank") %in% colnames(taxa))) {
    return(character())
  }
  query_names <- trimws(as.character(taxa$taxon_name))
  query_ranks <- as.character(taxa$taxon_rank)
  valid <- .is_non_empty_value(query_names) & .is_non_empty_value(query_ranks)
  if (!any(valid)) {
    return(character())
  }
  out <- query_ranks[valid]
  names(out) <- tolower(query_names[valid])
  out[!duplicated(names(out))]
}

.worms_query_rank <- function(query_key, query_rank_lookup) {
  if (length(query_rank_lookup) == 0) {
    return("unknown")
  }
  rank <- unname(query_rank_lookup[[tolower(as.character(query_key))]])
  if (length(rank) == 0 || is.null(rank) || !.is_non_empty_value(rank)) {
    return("unknown")
  }
  rank
}

.worms_query_key <- function(query, by) {
  if (identical(by, "aphia_id")) {
    return(.worms_aphia_id_key(query))
  }
  tolower(trimws(as.character(query)))
}

.validate_worms_taxa <- function(taxa, by) {
  if (missing(taxa) || is.null(taxa)) {
    stop("`taxa` must be supplied.", call. = FALSE)
  }
  if (inherits(taxa, "data.frame")) {
    if (!identical(by, "name")) {
      stop(
        "Data frame `taxa` input is only supported when `by = \"name\"`.",
        call. = FALSE
      )
    }
    if (!("taxon_name" %in% colnames(taxa))) {
      stop(
        "Data frame `taxa` input must contain a `taxon_name` column.",
        call. = FALSE
      )
    }
    taxa <- as.character(taxa$taxon_name)
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

  if (identical(by, "aphia_id")) {
    return(.normalise_worms_aphia_ids(
      taxa,
      drop_invalid = TRUE,
      warn = TRUE,
      value_name = "`taxa`"
    ))
  }

  taxa <- trimws(as.character(taxa))
  invalid <- is.na(taxa) | !nzchar(taxa)
  if (any(invalid)) {
    warning(
      "Dropped ",
      sum(invalid),
      " missing or empty taxon name value(s) from `taxa`.",
      call. = FALSE
    )
    taxa <- taxa[!invalid]
  }

  taxa
}

.deduplicate_worms_taxa <- function(taxa, by) {
  if (length(taxa) == 0) {
    return(taxa)
  }
  if (identical(by, "aphia_id")) {
    return(unique(taxa))
  }

  unique(taxa)
}

.normalise_worms_aphia_ids <- function(
  ids,
  drop_invalid = FALSE,
  warn = FALSE,
  value_name = "AphiaID values",
  empty_error = NULL
) {
  coerced <- .coerce_worms_aphia_ids(ids)
  invalid <- !coerced$valid
  if (any(invalid) && !isTRUE(drop_invalid)) {
    invalid_values <- unique(coerced$original[invalid])
    invalid_values[is.na(invalid_values)] <- "<NA>"
    stop(
      value_name,
      " must contain valid positive whole-number AphiaID values. Invalid ",
      "value(s): ",
      paste(invalid_values, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(invalid) && isTRUE(warn)) {
    warning(
      "Dropped ",
      sum(invalid),
      " missing or invalid AphiaID value(s) from ",
      value_name,
      ".",
      call. = FALSE
    )
  }

  out <- unique(unname(coerced$id[coerced$valid]))
  if (length(out) == 0) {
    if (is.null(empty_error)) {
      empty_error <- paste0(
        value_name,
        " must contain at least one valid positive whole-number AphiaID value."
      )
    }
    stop(empty_error, call. = FALSE)
  }

  out
}

.coerce_worms_aphia_ids <- function(ids) {
  original <- trimws(as.character(ids))
  original[is.na(ids)] <- NA_character_
  numeric_ids <- suppressWarnings(as.numeric(original))
  valid <- !is.na(original) &
    nzchar(original) &
    !is.na(numeric_ids) &
    is.finite(numeric_ids) &
    numeric_ids > 0 &
    numeric_ids == floor(numeric_ids)

  list(
    original = original,
    id = numeric_ids,
    valid = valid
  )
}

.worms_aphia_id_key <- function(ids) {
  coerced <- .coerce_worms_aphia_ids(ids)
  out <- rep(NA_character_, length(ids))
  out[coerced$valid] <- format(
    coerced$id[coerced$valid],
    scientific = FALSE,
    trim = TRUE
  )

  out
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
