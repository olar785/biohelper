#' Build a local WoRMS evidence database
#'
#' @description
#' Builds a reusable local WoRMS evidence database from a local WoRMS taxonomy
#' table or file. This is intended for once-off or occasional local database
#' builds, not for downloading the full WoRMS database and not for use inside
#' `flag_taxa()`.
#'
#' Users should obtain WoRMS taxonomy content according to the WoRMS terms of
#' use, then pass that local table to `build_worms_evidence_db()`. The function
#' enriches the supplied table by querying the WoRMS API through
#' `fetch_worms_evidence()` for taxonomic/environment evidence, saves the
#' standard evidence table as an RDS file, and can resume interrupted builds.
#'
#' All rows supplied by the user are considered. If you want to build evidence
#' for a rank or subset such as species only, filter the WoRMS taxonomy table
#' before calling this function. `build_worms_evidence_db()` does not interpret
#' evidence or make retain/review/exclude decisions.
#'
#' Large builds may take a long time. Use `batch_size`, `sleep`, and `resume =
#' TRUE` to query gently and cache progress.
#'
#' @param worms_taxonomy A data frame or path to a local WoRMS taxonomy file.
#'   Supported file types are `.csv`, `.tsv`, `.txt`, `.rds`, and `.xlsx`/`.xls`
#'   when the optional `readxl` package is installed.
#' @param output_path Path where the standard taxon evidence RDS file should be
#'   written.
#' @param id_col Optional AphiaID-like column name. Defaults to the first
#'   available of `AphiaID`, `aphia_id`, `aphiaID`, `taxonID`, `taxon_id`, or
#'   `source_taxon_id`. If an AphiaID-like column is available, AphiaID queries
#'   are always used. Missing or invalid AphiaID values are dropped before
#'   querying.
#' @param name_col Optional scientific-name column name. Defaults to the first
#'   available of `ScientificName`, `scientificName`, `scientificname`, `name`,
#'   or `taxon_name`. Names are used only when no AphiaID-like column is
#'   available.
#' @param batch_size Positive integer scalar giving the number of unique taxa to
#'   query per batch.
#' @param sleep Non-negative numeric scalar. Seconds to pause between batches
#'   and passed to `fetch_worms_evidence()` for AphiaID calls.
#' @param resume Logical scalar. If `TRUE` and `output_path` exists, read the
#'   existing evidence database and query only taxa not already present.
#' @param overwrite Logical scalar. If `TRUE`, rebuild from scratch even when
#'   `output_path` exists.
#' @param checked_at Date or character scalar recorded in new evidence rows.
#' @param verbose Logical scalar. If `TRUE`, emit progress messages.
#'
#' @return A standard taxon evidence data frame, also saved to `output_path` as
#'   an RDS file. The returned table is the combined full output from
#'   `fetch_worms_evidence()` across batches.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' evidence <- build_worms_evidence_db(
#'   worms_taxonomy = "worms_taxonomy_export.csv",
#'   output_path = "worms_evidence_db.rds",
#'   batch_size = 100,
#'   resume = TRUE
#' )
#' }
build_worms_evidence_db <- function(
  worms_taxonomy,
  output_path,
  id_col = NULL,
  name_col = NULL,
  batch_size = 100,
  sleep = 0.2,
  resume = TRUE,
  overwrite = FALSE,
  checked_at = Sys.Date(),
  verbose = TRUE
) {
  output_path <- .validate_required_path(output_path, "output_path")
  batch_size <- .validate_positive_integer_scalar(batch_size, "batch_size")
  sleep <- .validate_non_negative_numeric_scalar(sleep, "sleep")
  resume <- .validate_logical_scalar(resume, "resume")
  overwrite <- .validate_logical_scalar(overwrite, "overwrite")
  checked_at <- .validate_checked_at(checked_at)
  verbose <- .validate_logical_scalar(verbose, "verbose")

  if (file.exists(output_path) && !isTRUE(overwrite) && !isTRUE(resume)) {
    stop(
      "`output_path` already exists. Use `resume = TRUE` to continue from it, ",
      "or `overwrite = TRUE` to rebuild from scratch.",
      call. = FALSE
    )
  }

  taxonomy <- .read_worms_taxonomy_input(worms_taxonomy)
  query <- .extract_worms_taxonomy_queries(
    taxonomy = taxonomy,
    id_col = id_col,
    name_col = name_col
  )

  existing <- .initial_worms_evidence_db(
    output_path = output_path,
    resume = resume,
    overwrite = overwrite
  )
  query$cached_index <- .match_worms_cache(query$taxon, query$by[[1]], existing)
  taxa_to_query <- query$taxon[is.na(query$cached_index)]

  if (isTRUE(verbose)) {
    message(
      "Prepared ",
      nrow(query),
      " unique WoRMS taxon/taxa from local taxonomy input."
    )
    message(
      "Using ",
      if (identical(query$by[[1]], "aphia_id")) "AphiaID" else "scientific-name",
      " queries."
    )
    if (nrow(existing) > 0) {
      message(
        "Loaded ",
        nrow(existing),
        " existing WoRMS evidence row(s) from `output_path`."
      )
    }
    message(length(taxa_to_query), " taxon/taxa need new WoRMS queries.")
  }

  evidence <- existing
  if (length(taxa_to_query) == 0) {
    evidence <- validate_taxon_evidence(evidence)
    .write_worms_evidence_db(output_path, evidence)
    return(evidence)
  }

  batches <- split(taxa_to_query, ceiling(seq_along(taxa_to_query) / batch_size))
  for (batch_index in seq_along(batches)) {
    batch <- batches[[batch_index]]
    if (isTRUE(verbose)) {
      message(
        "Querying WoRMS batch ",
        batch_index,
        " of ",
        length(batches),
        " (",
        length(batch),
        " taxon/taxa)."
      )
    }

    batch_evidence <- fetch_worms_evidence(
      taxa = batch,
      by = query$by[[1]],
      marine_only = FALSE,
      cache_path = NULL,
      sleep = sleep,
      checked_at = checked_at,
      verbose = verbose
    )
    evidence <- .combine_taxon_evidence_rows(evidence, batch_evidence)
    .write_worms_evidence_db(output_path, evidence)

    if (batch_index < length(batches) && sleep > 0) {
      Sys.sleep(sleep)
    }
  }

  evidence
}

.read_worms_taxonomy_input <- function(worms_taxonomy) {
  if (inherits(worms_taxonomy, "data.frame")) {
    return(as.data.frame(
      worms_taxonomy,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  path <- .validate_required_path(worms_taxonomy, "worms_taxonomy")
  if (!file.exists(path)) {
    stop("`worms_taxonomy` file does not exist: ", path, call. = FALSE)
  }

  extension <- tolower(tools::file_ext(path))
  switch(
    extension,
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    txt = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    rds = {
      out <- readRDS(path)
      if (!inherits(out, "data.frame")) {
        stop("`worms_taxonomy` RDS file must contain a data.frame.", call. = FALSE)
      }
      as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    },
    xlsx = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop(
          "Reading Excel files requires the `readxl` package. ",
          "Install it with `install.packages(\"readxl\")`, or provide a ",
          ".csv, .tsv, .txt, or .rds file.",
          call. = FALSE
        )
      }
      as.data.frame(
        readxl::read_excel(path),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    },
    xls = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop(
          "Reading Excel files requires the `readxl` package. ",
          "Install it with `install.packages(\"readxl\")`, or provide a ",
          ".csv, .tsv, .txt, or .rds file.",
          call. = FALSE
        )
      }
      as.data.frame(
        readxl::read_excel(path),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    },
    stop(
      "`worms_taxonomy` must be a .csv, .tsv, .txt, .rds, .xlsx, or .xls file.",
      call. = FALSE
    )
  )
}

.extract_worms_taxonomy_queries <- function(
  taxonomy,
  id_col = NULL,
  name_col = NULL
) {
  if (!inherits(taxonomy, "data.frame")) {
    stop("`worms_taxonomy` must be a data.frame or supported file path.", call. = FALSE)
  }
  if (nrow(taxonomy) == 0 || ncol(taxonomy) == 0) {
    stop("`worms_taxonomy` must contain at least one row and one column.", call. = FALSE)
  }

  id_col <- .resolve_worms_taxonomy_column(
    taxonomy,
    id_col,
    c("AphiaID", "aphia_id", "aphiaID", "taxonID", "taxon_id", "source_taxon_id"),
    "id_col"
  )
  name_col <- .resolve_worms_taxonomy_column(
    taxonomy,
    name_col,
    c("ScientificName", "scientificName", "scientificname", "name", "taxon_name"),
    "name_col"
  )

  if (!is.null(id_col)) {
    query_col <- id_col
    by <- "aphia_id"
  } else if (!is.null(name_col)) {
    query_col <- name_col
    by <- "name"
  } else {
    stop(
      "`worms_taxonomy` must contain an AphiaID-like column or a scientific ",
      "name column.",
      call. = FALSE
    )
  }

  if (identical(by, "aphia_id")) {
    taxa <- .normalise_worms_aphia_ids(
      taxonomy[[query_col]],
      drop_invalid = TRUE,
      warn = TRUE,
      value_name = paste0("`", query_col, "`"),
      empty_error = paste0(
        "All AphiaID values in `",
        query_col,
        "` are missing or invalid after coercion. Check `id_col`, or remove ",
        "the AphiaID-like column to use a scientific-name fallback."
      )
    )
  } else {
    taxa <- trimws(as.character(taxonomy[[query_col]]))
    taxa <- taxa[.is_non_empty_value(taxa)]
    taxa <- unique(taxa)
    if (length(taxa) == 0) {
      stop("No non-empty WoRMS taxa were found in `worms_taxonomy`.", call. = FALSE)
    }
  }

  data.frame(
    taxon = taxa,
    by = by,
    query_col = query_col,
    id_col = if (is.null(id_col)) NA_character_ else id_col,
    name_col = if (is.null(name_col)) NA_character_ else name_col,
    stringsAsFactors = FALSE
  )
}

.resolve_worms_taxonomy_column <- function(taxonomy, column, candidates, argument_name) {
  if (!is.null(column)) {
    column <- .validate_required_scalar_character(column, argument_name)
    matches <- which(tolower(colnames(taxonomy)) == tolower(column))
    if (length(matches) == 0) {
      stop(
        "`",
        argument_name,
        "` was not found in `worms_taxonomy`: ",
        column,
        call. = FALSE
      )
    }
    return(colnames(taxonomy)[[matches[[1]]]])
  }

  matches <- match(tolower(candidates), tolower(colnames(taxonomy)))
  matches <- matches[!is.na(matches)]
  if (length(matches) == 0) {
    return(NULL)
  }

  colnames(taxonomy)[[matches[[1]]]]
}

.initial_worms_evidence_db <- function(output_path, resume, overwrite) {
  if (isTRUE(overwrite) || !file.exists(output_path)) {
    return(.empty_taxon_evidence_table())
  }
  if (!isTRUE(resume)) {
    return(.empty_taxon_evidence_table())
  }

  evidence <- readRDS(output_path)
  if (!inherits(evidence, "data.frame")) {
    stop("Existing `output_path` must contain a data.frame saved as RDS.", call. = FALSE)
  }

  validate_taxon_evidence(evidence)
}

.write_worms_evidence_db <- function(output_path, evidence) {
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  evidence <- validate_taxon_evidence(evidence)
  saveRDS(evidence, output_path)
  invisible(evidence)
}

.combine_taxon_evidence_rows <- function(...) {
  evidence_list <- list(...)
  evidence_list <- Filter(Negate(is.null), evidence_list)
  if (length(evidence_list) == 0) {
    return(.empty_taxon_evidence_table())
  }

  evidence_list <- lapply(evidence_list, function(evidence) {
    if (!inherits(evidence, "data.frame")) {
      stop("WoRMS evidence batches must be data.frames.", call. = FALSE)
    }
    validate_taxon_evidence(evidence)
  })

  all_columns <- unique(unlist(lapply(evidence_list, colnames), use.names = FALSE))
  evidence_list <- lapply(evidence_list, function(evidence) {
    missing_columns <- setdiff(all_columns, colnames(evidence))
    if (length(missing_columns) > 0) {
      evidence[missing_columns] <- rep(
        list(rep(NA_character_, nrow(evidence))),
        length(missing_columns)
      )
    }
    evidence[, all_columns, drop = FALSE]
  })

  out <- do.call(rbind, evidence_list)
  rownames(out) <- NULL
  validate_taxon_evidence(out)
}

.validate_required_path <- function(value, name) {
  if (
    missing(value) ||
      !is.character(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !nzchar(trimws(value))
  ) {
    stop("`", name, "` must be a non-empty character scalar.", call. = FALSE)
  }

  trimws(value)
}

.validate_positive_integer_scalar <- function(value, name) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      value < 1 ||
      value != as.integer(value)
  ) {
    stop("`", name, "` must be a positive integer scalar.", call. = FALSE)
  }

  as.integer(value)
}
