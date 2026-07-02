#' Build a standard taxon evidence table
#'
#' @description
#' Creates a local, standardised taxon evidence table from user-provided source
#' tables or database exports. The result can be passed to
#' `flag_taxa(taxon_evidence = ...)` to provide reusable ecological,
#' geographic, occurrence, barcode, or taxonomic database evidence.
#'
#' `build_taxon_evidence()` does not download full source databases and does not
#' make live API calls. Users should obtain source tables externally, for
#' example from WoRMS, OBIS, GBIF, BOLD, or curated literature summaries, and
#' pass those tables via `sources` or `user_evidence`.
#'
#' @param sources A named list of data frames, such as
#'   `list(worms = worms_df, obis = obis_df, gbif = gbif_df, bold = bold_df)`.
#'   A single data frame may be supplied when `source_type` is provided.
#' @param user_evidence Optional user-supplied evidence table already in the
#'   standard format accepted by `flag_taxa(taxon_evidence = ...)`.
#' @param source_type Optional character vector describing how to interpret an
#'   unnamed source input. Supported values are `"worms"`, `"obis"`, `"gbif"`,
#'   and `"bold"`.
#' @param verbose Logical scalar. If `TRUE`, emit a short message about the
#'   number of evidence rows returned.
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
#' `basis_of_record`, `occurrence_count`, `reference_url`, `doi`,
#' `checked_at`, and lineage columns such as `kingdom`, `phylum`, `class`,
#' `order`, `family`, and `genus`.
#'
#' @section Supported source tables:
#' Source-specific import helpers map common export columns when present and
#' fill missing standard columns with `NA`. These helpers intentionally support
#' only simple, common source-table shapes at this stage.
#'
#' WoRMS-like tables can include columns such as `scientificname`,
#' `scientificName`, `valid_name`, `accepted_name`, `rank`, `taxonRank`,
#' `AphiaID`, `valid_AphiaID`, `isMarine`, `isBrackish`, `isFreshwater`,
#' `isTerrestrial`, `url`, or `lsid`.
#'
#' OBIS-like and GBIF-like tables can include occurrence columns such as
#' `scientificName`, `taxonRank`, source taxon identifiers, occurrence
#' identifiers, coordinates, locality, country, `basisOfRecord`, and
#' `individualCount`.
#'
#' BOLD-like tables can include taxonomic columns such as `species_name`,
#' `genus_name`, `family_name`, `order_name`, `class_name`, and `phylum_name`,
#' plus specimen fields such as `processid`, `sampleid`, `country`,
#' `province_state`, `region`, `exactsite`, `lat`, and `lon`.
#'
#' @export
#'
#' @examples
#' worms_export <- data.frame(
#'   scientificname = "Salmo salar",
#'   rank = "Species",
#'   AphiaID = "127186",
#'   isMarine = TRUE
#' )
#'
#' evidence <- build_taxon_evidence(sources = list(worms = worms_export))
build_taxon_evidence <- function(
  sources = list(),
  user_evidence = NULL,
  source_type = NULL,
  verbose = FALSE
) {
  verbose <- .validate_logical_scalar(verbose, "verbose")
  sources <- .validate_taxon_evidence_sources(sources, source_type)
  user_evidence <- validate_taxon_evidence(user_evidence)

  source_evidence <- lapply(names(sources), function(source_name) {
    .standardise_source_evidence(sources[[source_name]], source_name)
  })

  evidence_parts <- c(
    lapply(source_evidence, .standardise_taxon_evidence_columns),
    list(.standardise_taxon_evidence_columns(user_evidence))
  )
  evidence <- do.call(rbind, evidence_parts)
  evidence <- validate_taxon_evidence(evidence)

  if (isTRUE(verbose)) {
    message("Built ", nrow(evidence), " taxon evidence row(s).")
  }

  evidence
}

standardise_worms_evidence <- function(x) {
  x <- .validate_source_data_frame(x, "x")
  out <- .empty_taxon_evidence_table(nrow(x))
  if (nrow(x) == 0) {
    return(out)
  }

  source_taxon_id <- .row_first_non_empty(x, c(
    "valid_AphiaID",
    "AphiaID",
    "aphia_id"
  ))
  environment <- .worms_environment_summary(x)

  out$taxon_name <- .row_first_non_empty(x, c(
    "scientificname",
    "scientificName",
    "valid_name",
    "accepted_name"
  ))
  out$taxon_rank <- .row_first_non_empty(x, c("rank", "taxonRank"))
  out$source <- "worms"
  out$evidence_type <- "taxonomic_database"
  out$evidence_summary <- .database_evidence_summary(
    "WoRMS taxonomic database record",
    environment
  )
  out$reference <- "WoRMS"
  out$accepted_name <- .row_first_non_empty(x, c(
    "valid_name",
    "accepted_name",
    "scientificname",
    "scientificName"
  ))
  out$accepted_rank <- .row_first_non_empty(x, c("rank", "taxonRank"))
  out$source_taxon_id <- source_taxon_id
  out$environment <- environment
  out$reference_url <- .worms_reference_url(x, source_taxon_id)
  out$checked_at <- as.character(Sys.Date())

  out
}

standardise_obis_evidence <- function(x) {
  x <- .validate_source_data_frame(x, "x")
  out <- .empty_taxon_evidence_table(nrow(x))
  if (nrow(x) == 0) {
    return(out)
  }

  out$taxon_name <- .row_first_non_empty(x, c(
    "scientificName",
    "scientificname",
    "acceptedScientificName"
  ))
  out$taxon_rank <- .row_first_non_empty(x, c("taxonRank", "rank"))
  out$source <- "obis"
  out$evidence_type <- "occurrence"
  out$evidence_summary <- "OBIS occurrence record."
  out$reference <- "OBIS"
  out$accepted_name <- .row_first_non_empty(x, c(
    "acceptedScientificName",
    "scientificName",
    "scientificname"
  ))
  out$accepted_rank <- out$taxon_rank
  out$source_taxon_id <- .row_first_non_empty(x, c("AphiaID", "aphia_id"))
  out$source_record_id <- .row_first_non_empty(x, c("occurrenceID", "id"))
  out$environment <- "marine"
  out$region <- .row_first_non_empty(x, c("country", "countryCode", "region"))
  out$locality <- .row_first_non_empty(x, c(
    "locality",
    "area",
    "areaName",
    "locationID"
  ))
  out$decimal_latitude <- .row_first_non_empty(x, c(
    "decimalLatitude",
    "decimal_latitude"
  ))
  out$decimal_longitude <- .row_first_non_empty(x, c(
    "decimalLongitude",
    "decimal_longitude"
  ))
  out$basis_of_record <- .row_first_non_empty(x, c("basisOfRecord", "basis_of_record"))
  out$occurrence_count <- .row_first_non_empty(x, c(
    "individualCount",
    "occurrenceCount",
    "occurrence_count"
  ))
  out$reference_url <- .row_first_non_empty(x, c("references", "url", "reference_url"))
  out$doi <- .row_first_non_empty(x, c("doi", "DOI"))
  out$checked_at <- as.character(Sys.Date())

  out
}

standardise_gbif_evidence <- function(x) {
  x <- .validate_source_data_frame(x, "x")
  out <- .empty_taxon_evidence_table(nrow(x))
  if (nrow(x) == 0) {
    return(out)
  }

  out$taxon_name <- .row_first_non_empty(x, c(
    "scientificName",
    "scientificname",
    "acceptedScientificName",
    "canonicalName"
  ))
  out$taxon_rank <- .row_first_non_empty(x, c("taxonRank", "rank"))
  out$source <- "gbif"
  out$evidence_type <- "occurrence"
  out$evidence_summary <- "GBIF occurrence record."
  out$reference <- "GBIF"
  out$accepted_name <- .row_first_non_empty(x, c(
    "acceptedScientificName",
    "scientificName",
    "scientificname",
    "canonicalName"
  ))
  out$accepted_rank <- out$taxon_rank
  out$source_taxon_id <- .row_first_non_empty(x, c("taxonKey", "usageKey"))
  out$source_record_id <- .row_first_non_empty(x, c("gbifID", "occurrenceID"))
  out$region <- .row_first_non_empty(x, c("country", "countryCode", "region"))
  out$locality <- .row_first_non_empty(x, c("locality", "stateProvince", "county"))
  out$decimal_latitude <- .row_first_non_empty(x, c(
    "decimalLatitude",
    "decimal_latitude"
  ))
  out$decimal_longitude <- .row_first_non_empty(x, c(
    "decimalLongitude",
    "decimal_longitude"
  ))
  out$basis_of_record <- .row_first_non_empty(x, c("basisOfRecord", "basis_of_record"))
  out$occurrence_count <- .row_first_non_empty(x, c("individualCount", "occurrenceCount"))
  out$reference_url <- .row_first_non_empty(x, c("references", "url", "reference_url"))
  out$doi <- .row_first_non_empty(x, c("doi", "DOI"))
  out$checked_at <- as.character(Sys.Date())

  out
}

standardise_bold_evidence <- function(x) {
  x <- .validate_source_data_frame(x, "x")
  out <- .empty_taxon_evidence_table(nrow(x))
  if (nrow(x) == 0) {
    return(out)
  }

  bold_taxa <- .bold_query_taxa(x)

  out$taxon_name <- bold_taxa$taxon_name
  out$taxon_rank <- bold_taxa$taxon_rank
  out$source <- "bold"
  out$evidence_type <- "barcode_or_specimen"
  out$evidence_summary <- "BOLD barcode or specimen record."
  out$reference <- "BOLD"
  out$accepted_name <- bold_taxa$taxon_name
  out$accepted_rank <- bold_taxa$taxon_rank
  out$source_record_id <- .row_first_non_empty(x, c("processid", "sampleid"))
  out$region <- .row_join_non_empty(x, c("country", "province_state", "region"))
  out$locality <- .row_first_non_empty(x, c("exactsite", "locality"))
  out$decimal_latitude <- .row_first_non_empty(x, c("lat", "latitude"))
  out$decimal_longitude <- .row_first_non_empty(x, c("lon", "lng", "longitude"))
  out$checked_at <- as.character(Sys.Date())

  out
}

.validate_taxon_evidence_sources <- function(sources, source_type = NULL) {
  if (is.null(sources)) {
    sources <- list()
  }

  if (inherits(sources, "data.frame")) {
    if (is.null(source_type)) {
      stop(
        "`source_type` must be provided when `sources` is a single data.frame.",
        call. = FALSE
      )
    }
    source_type <- .validate_character_vector(source_type, "source_type")
    if (length(source_type) != 1) {
      stop(
        "`source_type` must contain exactly one value when `sources` is a ",
        "single data.frame.",
        call. = FALSE
      )
    }
    sources <- stats::setNames(list(sources), source_type)
  }

  if (!is.list(sources)) {
    stop("`sources` must be a named list of data.frames.", call. = FALSE)
  }
  if (length(sources) == 0) {
    return(sources)
  }

  source_type <- .validate_optional_source_type(source_type)
  source_names <- names(sources)
  if (is.null(source_names)) {
    source_names <- rep("", length(sources))
  }

  if (length(sources) == 1 && !is.null(source_type)) {
    source_names <- source_type[[1]]
  } else {
    source_names <- .fill_missing_source_names(source_names, source_type)
  }

  source_names <- tolower(trimws(source_names))
  unsupported_sources <- setdiff(unique(source_names), .supported_evidence_sources())
  if (length(unsupported_sources) > 0) {
    stop(
      "Unsupported `sources` names: ",
      paste(unsupported_sources, collapse = ", "),
      ". Supported values are: ",
      paste(.supported_evidence_sources(), collapse = ", "),
      call. = FALSE
    )
  }

  for (source_index in seq_along(sources)) {
    source_name <- source_names[[source_index]]
    .validate_source_data_frame(
      sources[[source_index]],
      paste0("sources$", source_name)
    )
  }

  names(sources) <- source_names
  sources
}

.validate_optional_source_type <- function(source_type) {
  if (is.null(source_type)) {
    return(NULL)
  }

  .validate_character_vector(source_type, "source_type")
}

.fill_missing_source_names <- function(source_names, source_type) {
  missing_names <- !nzchar(trimws(source_names))
  if (!any(missing_names)) {
    return(source_names)
  }
  if (is.null(source_type)) {
    stop(
      "`sources` must be a named list of data.frames, or `source_type` must ",
      "be provided for unnamed source inputs.",
      call. = FALSE
    )
  }

  if (length(source_type) == sum(missing_names)) {
    source_names[missing_names] <- source_type
    return(source_names)
  }
  if (length(source_type) == length(source_names)) {
    source_names[missing_names] <- source_type[missing_names]
    return(source_names)
  }

  stop(
    "`source_type` must have length 1 for a single unnamed source, length ",
    "matching the unnamed sources, or length matching all sources.",
    call. = FALSE
  )
}

.validate_source_data_frame <- function(x, name) {
  if (!inherits(x, "data.frame")) {
    stop("`", name, "` must be a data.frame.", call. = FALSE)
  }

  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

.supported_evidence_sources <- function() {
  c("worms", "obis", "gbif", "bold")
}

.standardise_source_evidence <- function(x, source_name) {
  switch(
    source_name,
    worms = standardise_worms_evidence(x),
    obis = standardise_obis_evidence(x),
    gbif = standardise_gbif_evidence(x),
    bold = standardise_bold_evidence(x)
  )
}

.source_column <- function(x, candidates) {
  column_index <- match(tolower(candidates), tolower(colnames(x)))
  column_index <- column_index[!is.na(column_index)]
  if (length(column_index) == 0) {
    return(rep(NA_character_, nrow(x)))
  }

  as.character(x[[column_index[[1]]]])
}

.has_source_column <- function(x, candidates) {
  any(tolower(candidates) %in% tolower(colnames(x)))
}

.row_first_non_empty <- function(x, candidates) {
  columns <- lapply(candidates, function(candidate) .source_column(x, candidate))
  values <- do.call(cbind, columns)
  out <- apply(values, 1, function(row_values) {
    row_values <- row_values[.is_non_empty_value(row_values)]
    if (length(row_values) == 0) {
      return(NA_character_)
    }
    row_values[[1]]
  })
  unname(as.character(out))
}

.row_join_non_empty <- function(x, candidates, collapse = "; ") {
  columns <- lapply(candidates, function(candidate) .source_column(x, candidate))
  values <- do.call(cbind, columns)
  out <- apply(values, 1, function(row_values) {
    row_values <- unique(row_values[.is_non_empty_value(row_values)])
    if (length(row_values) == 0) {
      return(NA_character_)
    }
    paste(row_values, collapse = collapse)
  })
  unname(as.character(out))
}

.is_non_empty_value <- function(value) {
  !is.na(value) & nzchar(trimws(as.character(value)))
}

.truthy_source_value <- function(value) {
  tolower(trimws(as.character(value))) %in% c("1", "true", "t", "yes", "y")
}

.worms_environment_summary <- function(x) {
  environment_columns <- c(
    isMarine = "marine",
    isBrackish = "brackish",
    isFreshwater = "freshwater",
    isTerrestrial = "terrestrial"
  )
  available_columns <- names(environment_columns)[
    vapply(
      names(environment_columns),
      function(column_name) .has_source_column(x, column_name),
      logical(1)
    )
  ]

  if (length(available_columns) == 0) {
    return(rep(NA_character_, nrow(x)))
  }

  flags <- lapply(available_columns, function(column_name) {
    .truthy_source_value(.source_column(x, column_name))
  })
  names(flags) <- unname(environment_columns[available_columns])

  out <- vapply(seq_len(nrow(x)), function(row_index) {
    environments <- names(flags)[
      vapply(flags, function(flag) isTRUE(flag[[row_index]]), logical(1))
    ]
    if (length(environments) == 0) {
      return(NA_character_)
    }
    paste(environments, collapse = "; ")
  }, character(1))

  out
}

.database_evidence_summary <- function(prefix, environment) {
  out <- paste0(prefix, ".")
  has_environment <- .is_non_empty_value(environment)
  out[has_environment] <- paste0(
    prefix,
    "; environment flags: ",
    environment[has_environment],
    "."
  )

  out
}

.worms_reference_url <- function(x, source_taxon_id) {
  reference_url <- .row_first_non_empty(x, c("url", "URL", "lsid", "LSID"))
  missing_url <- !.is_non_empty_value(reference_url) &
    .is_non_empty_value(source_taxon_id)
  reference_url[missing_url] <- paste0(
    "https://www.marinespecies.org/aphia.php?p=taxdetails&id=",
    source_taxon_id[missing_url]
  )

  reference_url
}

.bold_query_taxa <- function(x) {
  rank_columns <- c(
    species = "species_name",
    genus = "genus_name",
    family = "family_name",
    order = "order_name",
    class = "class_name",
    phylum = "phylum_name"
  )
  taxon_name <- rep(NA_character_, nrow(x))
  taxon_rank <- rep(NA_character_, nrow(x))

  for (rank in names(rank_columns)) {
    values <- .source_column(x, rank_columns[[rank]])
    use_value <- !.is_non_empty_value(taxon_name) & .is_non_empty_value(values)
    taxon_name[use_value] <- values[use_value]
    taxon_rank[use_value] <- rank
  }

  data.frame(
    taxon_name = taxon_name,
    taxon_rank = taxon_rank,
    stringsAsFactors = FALSE
  )
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
