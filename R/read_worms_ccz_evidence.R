#' Read local WoRMS Clarion-Clipperton Zone checklist evidence
#'
#' @description
#' Reads a locally downloaded WoRMS Clarion-Clipperton Zone (CCZ)
#' taxlist/checklist file and converts it to the standard `taxon_evidence`
#' format used by `flag_taxa()`. This function is a local file reader and
#' standardiser only: it does not download data, query WoRMS, or call any LLM.
#'
#' Users must obtain the CCZ checklist separately from WoRMS according to the
#' WoRMS terms of use. The returned rows use `source = "worms_ccz"` and
#' `evidence_type = "regional_deepsea_checklist"` so this regional deep-sea
#' evidence remains separate from generic WoRMS evidence returned by
#' `fetch_worms_evidence()`. Both evidence types can be combined and passed to
#' `flag_taxa(taxon_evidence = ...)`.
#'
#' @param path Path to a local CCZ checklist file. Supported file types are
#'   `.csv`, `.tsv`, `.txt`, `.rds`, and `.xlsx`/`.xls` when the optional
#'   `readxl` package is installed.
#' @param checked_at Date or character scalar recorded in the `checked_at`
#'   column.
#'
#' @return A data frame in the standard taxon evidence format, with optional
#'   lineage columns (`kingdom`, `phylum`, `class`, `order`, `family`, and
#'   `genus`) preserved when present in the input.
#' @export
#'
#' @examples
#' \dontrun{
#' ccz_evidence <- read_worms_ccz_evidence("CCZ_taxlist.csv")
#'
#' taxa_to_query <- extract_taxa_for_evidence(ps_test_data_euk)
#' worms_evidence <- fetch_worms_evidence(
#'   taxa_to_query,
#'   by = "name",
#'   cache_path = "worms_cache.rds"
#' )
#'
#' taxon_evidence <- dplyr::bind_rows(worms_evidence, ccz_evidence)
#' prompt <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   expected_habitat = "deep sea",
#'   expected_region = "Clarion-Clipperton Zone",
#'   taxon_evidence = taxon_evidence,
#'   prompt_only = TRUE
#' )
#' }
read_worms_ccz_evidence <- function(path, checked_at = Sys.Date()) {
  path <- .validate_required_path(path, "path")
  checked_at <- .validate_checked_at(checked_at)

  ccz_taxlist <- .read_worms_ccz_input(path)
  evidence <- .standardise_worms_ccz_evidence(ccz_taxlist, checked_at = checked_at)
  validate_taxon_evidence(evidence)
}

.read_worms_ccz_input <- function(path) {
  if (!file.exists(path)) {
    stop("`path` file does not exist: ", path, call. = FALSE)
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
        stop("`path` RDS file must contain a data.frame.", call. = FALSE)
      }
      as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    },
    xlsx = .read_worms_ccz_excel(path),
    xls = .read_worms_ccz_excel(path),
    stop(
      "`path` must be a .csv, .tsv, .txt, .rds, .xlsx, or .xls file.",
      call. = FALSE
    )
  )
}

.read_worms_ccz_excel <- function(path) {
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
}

.standardise_worms_ccz_evidence <- function(x, checked_at) {
  x <- .validate_source_data_frame(x, "path")
  out <- .empty_taxon_evidence_table(nrow(x))
  if (nrow(x) == 0) {
    return(out)
  }

  source_taxon_id <- .row_first_non_empty(x, c("AphiaID", "aphia_id", "source_taxon_id"))
  source_record_id <- .row_first_non_empty(x, c(
    "AphiaID_accepted",
    "valid_AphiaID",
    "source_record_id"
  ))
  taxon_rank <- .row_first_non_empty(x, c("taxonRank", "rank", "taxon_rank"))
  status <- .row_first_non_empty(x, c("taxonomicStatus", "status"))
  environment <- .worms_ccz_environment_summary(x)
  reference <- .row_first_non_empty(x, c("Citation", "citation", "reference"))
  missing_reference <- !.is_non_empty_value(reference)
  reference[missing_reference] <- "WoRMS Clarion-Clipperton Zone checklist"

  out$taxon_name <- .row_first_non_empty(x, c(
    "ScientificName",
    "scientificName",
    "scientificname",
    "taxon_name"
  ))
  out$taxon_rank <- taxon_rank
  out$source <- "worms_ccz"
  out$evidence_type <- "regional_deepsea_checklist"
  out$evidence_summary <- .worms_ccz_evidence_summary(status, environment)
  out$reference <- reference
  out$accepted_name <- .row_first_non_empty(x, c(
    "ScientificName_accepted",
    "scientificName_accepted",
    "valid_name",
    "accepted_name"
  ))
  out$accepted_rank <- taxon_rank
  out$source_taxon_id <- source_taxon_id
  out$source_record_id <- source_record_id
  out$environment <- environment
  out$habitat <- "deep sea"
  out$region <- "Clarion-Clipperton Zone"
  out$locality <- .row_first_non_empty(x, c("locality", "Locality"))
  out$decimal_latitude <- .row_first_non_empty(x, c(
    "decimalLatitude",
    "decimal_latitude"
  ))
  out$decimal_longitude <- .row_first_non_empty(x, c(
    "decimalLongitude",
    "decimal_longitude"
  ))
  out$basis_of_record <- .row_first_non_empty(x, c(
    "basisOfRecord",
    "basis_of_record"
  ))
  out$occurrence_count <- .row_first_non_empty(x, c(
    "occurrenceCount",
    "occurrence_count",
    "individualCount"
  ))
  out$reference_url <- .worms_ccz_reference_url(
    x = x,
    reference = reference,
    source_taxon_id = source_taxon_id
  )
  out$doi <- .row_first_non_empty(x, c("doi", "DOI"))
  out$checked_at <- checked_at
  out$kingdom <- .row_first_non_empty(x, c("Kingdom", "kingdom"))
  out$phylum <- .row_first_non_empty(x, c("Phylum", "phylum"))
  out$class <- .row_first_non_empty(x, c("Class", "class"))
  out$order <- .row_first_non_empty(x, c("Order", "order"))
  out$family <- .row_first_non_empty(x, c("Family", "family"))
  out$genus <- .row_first_non_empty(x, c("Genus", "genus"))

  validate_taxon_evidence(out)
}

.worms_ccz_environment_summary <- function(x) {
  environment_columns <- c(
    Marine = "marine",
    Brackish = "brackish",
    Fresh = "freshwater",
    Freshwater = "freshwater",
    Terrestrial = "terrestrial"
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

  vapply(seq_len(nrow(x)), function(row_index) {
    environments <- unique(names(flags)[
      vapply(flags, function(flag) isTRUE(flag[[row_index]]), logical(1))
    ])
    if (length(environments) == 0) {
      return(NA_character_)
    }
    paste(environments, collapse = "; ")
  }, character(1))
}

.worms_ccz_evidence_summary <- function(status, environment) {
  out <- rep(
    "Taxon is listed in the WoRMS Clarion-Clipperton Zone checklist.",
    length(status)
  )

  has_status <- .is_non_empty_value(status)
  out[has_status] <- paste0(
    out[has_status],
    " Taxonomic status: ",
    status[has_status],
    "."
  )

  has_environment <- .is_non_empty_value(environment)
  out[has_environment] <- paste0(
    out[has_environment],
    " Environment flags: ",
    environment[has_environment],
    "."
  )

  out
}

.worms_ccz_reference_url <- function(x, reference, source_taxon_id) {
  reference_url <- .extract_first_url(reference)
  source_url <- .row_first_non_empty(x, c("reference_url", "url", "URL"))
  missing_url <- !.is_non_empty_value(reference_url) & .is_non_empty_value(source_url)
  reference_url[missing_url] <- source_url[missing_url]

  missing_url <- !.is_non_empty_value(reference_url) & .is_non_empty_value(source_taxon_id)
  reference_url[missing_url] <- paste0(
    "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=",
    source_taxon_id[missing_url]
  )

  reference_url
}

.extract_first_url <- function(value) {
  value <- as.character(value)
  out <- rep(NA_character_, length(value))
  url_match <- regexpr("https?://[^[:space:]\"'<>]+", value)
  url_length <- attr(url_match, "match.length")
  has_url <- !is.na(url_match) & url_match > 0
  if (any(has_url)) {
    out[has_url] <- substring(
      value[has_url],
      url_match[has_url],
      url_match[has_url] + url_length[has_url] - 1
    )
  }
  out
}
