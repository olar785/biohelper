#' Flag taxa against expected ecological and geographic context
#'
#' @description
#' Builds a review prompt for checking whether taxa are compatible with an
#' expected environment, habitat, and region. This first scaffold does not call
#' an LLM or query external databases. It only extracts taxonomy, chooses the
#' most specific useful taxon for each row, and returns the prompt that can be
#' inspected or copied into a later review workflow. It can also validate a
#' mocked LLM-like result table while the backend is under development.
#'
#' @param x A taxonomy table as a data frame, matrix, phyloseq taxonomy table,
#'   or a phyloseq/speedyseq object containing a taxonomy table.
#' @param expected_environment Character scalar describing the expected broad
#'   environment, for example `"marine"`, `"freshwater"`, or `"terrestrial"`.
#' @param expected_habitat Optional character scalar describing the expected
#'   habitat, for example `"estuary"` or `"kelp forest"`.
#' @param expected_region Optional character scalar describing the expected
#'   geographic region.
#' @param tax_ranks Character vector of taxonomic ranks to inspect, ordered from
#'   most specific to broadest.
#' @param evidence_sources Character vector naming evidence sources the prompt
#'   should ask the reviewer to consider.
#' @param taxon_evidence Optional data frame containing user-supplied evidence
#'   about taxon environment, habitat, and/or region. Required columns are
#'   `taxon_name`, `taxon_rank`, `source`, `evidence_type`, `evidence_summary`,
#'   and `reference`. Optional columns are `accepted_name`, `accepted_rank`,
#'   `source_taxon_id`, `source_record_id`, `environment`, `habitat`, `region`,
#'   `locality`, `decimal_latitude`, `decimal_longitude`, `basis_of_record`,
#'   `occurrence_count`, `reference_url`, `doi`, and `checked_at`.
#' @param prompt_only Logical scalar. When `TRUE`, the generated prompt is
#'   returned as a character scalar. When `FALSE`, `mock_llm_result` must be
#'   provided because the LLM backend is not implemented yet.
#' @param verbose Logical scalar. If `TRUE`, emit a short message about the
#'   number of taxa prepared for the prompt.
#' @param mock_llm_result Optional data frame used for testing the validation
#'   of an LLM-like result while the backend is not implemented.
#'
#' @return A character scalar containing the generated prompt when
#'   `prompt_only = TRUE`, or a validated result table when `prompt_only = FALSE`
#'   and `mock_llm_result` is provided.
#'
#' @section Taxon evidence table:
#' `taxon_evidence` can be used to provide curated evidence that should be
#' considered before any future online evidence lookup. The required columns are:
#' `taxon_name`, `taxon_rank`, `source`, `evidence_type`, `evidence_summary`,
#' and `reference`.
#'
#' Optional columns are: `accepted_name`, `accepted_rank`, `source_taxon_id`,
#' `source_record_id`, `environment`, `habitat`, `region`, `locality`,
#' `decimal_latitude`, `decimal_longitude`, `basis_of_record`,
#' `occurrence_count`, `reference_url`, `doi`, and `checked_at`.
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
#' prompt <- flag_taxa(
#'   tax,
#'   expected_environment = "marine",
#'   expected_habitat = "coastal water",
#'   expected_region = "North Atlantic"
#' )
flag_taxa <- function(
  x,
  expected_environment,
  expected_habitat = NULL,
  expected_region = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum", "kingdom"),
  evidence_sources = c("worms", "obis", "gbif", "bold", "literature"),
  taxon_evidence = NULL,
  prompt_only = TRUE,
  verbose = FALSE,
  mock_llm_result = NULL
) {
  expected_environment <- .validate_required_scalar_character(
    expected_environment,
    "expected_environment"
  )
  expected_habitat <- .validate_optional_scalar_character(
    expected_habitat,
    "expected_habitat"
  )
  expected_region <- .validate_optional_scalar_character(
    expected_region,
    "expected_region"
  )
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")
  evidence_sources <- .validate_character_vector(
    evidence_sources,
    "evidence_sources"
  )
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  prompt_only <- .validate_logical_scalar(prompt_only, "prompt_only")
  verbose <- .validate_logical_scalar(verbose, "verbose")

  tax_table <- extract_tax_table(x)
  if (isTRUE(verbose)) {
    message("Prepared ", nrow(tax_table), " taxonomic rows for prompt-only review.")
  }

  if (!isTRUE(prompt_only) && !is.null(mock_llm_result)) {
    return(validate_flag_taxa_output(
      result = mock_llm_result,
      original_taxonomy = tax_table,
      expected_region = expected_region
    ))
  }

  if (!isTRUE(prompt_only)) {
    stop(
      "The LLM backend for `flag_taxa()` is not implemented yet. ",
      "Use `prompt_only = TRUE` to return the generated prompt, or provide ",
      "`mock_llm_result` for validation testing.",
      call. = FALSE
    )
  }

  build_flag_taxa_prompt(
    tax_table = tax_table,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    tax_ranks = tax_ranks,
    evidence_sources = evidence_sources,
    taxon_evidence = taxon_evidence
  )
}

#' Validate a flag_taxa result table
#'
#' @param result Data frame returned by an LLM-like review.
#' @param original_taxonomy Original taxonomy table used to build the prompt.
#' @param expected_region Optional expected region supplied to `flag_taxa()`.
#'
#' @return `result`, unchanged.
#' @noRd
validate_flag_taxa_output <- function(result, original_taxonomy, expected_region = NULL) {
  if (!inherits(result, "data.frame")) {
    stop("`result` must be a data.frame.", call. = FALSE)
  }
  if (!inherits(original_taxonomy, "data.frame")) {
    stop("`original_taxonomy` must be a data.frame.", call. = FALSE)
  }

  required_columns <- .flag_taxa_output_columns()
  missing_columns <- setdiff(required_columns, colnames(result))
  if (length(missing_columns) > 0) {
    stop(
      "`result` is missing required flag_taxa columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(result) != nrow(original_taxonomy)) {
    stop(
      "`result` must have the same number of rows as `original_taxonomy`. ",
      "`result` has ",
      nrow(result),
      " row(s); `original_taxonomy` has ",
      nrow(original_taxonomy),
      " row(s).",
      call. = FALSE
    )
  }

  if (
    "feature_id" %in% colnames(result) &&
      "feature_id" %in% colnames(original_taxonomy)
  ) {
    result_feature_id <- as.character(result$feature_id)
    original_feature_id <- as.character(original_taxonomy$feature_id)
    if (!identical(result_feature_id, original_feature_id)) {
      stop(
        "`feature_id` values in `result` must match `original_taxonomy` ",
        "in the same row order.",
        call. = FALSE
      )
    }
  }

  allowed_values <- .flag_taxa_allowed_values()
  for (column_name in names(allowed_values)) {
    .validate_flag_taxa_enum_column(
      result = result,
      column_name = column_name,
      allowed_values = allowed_values[[column_name]]
    )
  }
  .validate_flag_taxa_region_assessment(result, expected_region)

  result
}

#' Validate user-provided taxon evidence
#'
#' @param taxon_evidence Optional evidence table supplied to `flag_taxa()`.
#'
#' @return `NULL`, or a data frame with standardised column order.
#' @noRd
validate_taxon_evidence <- function(taxon_evidence) {
  if (is.null(taxon_evidence)) {
    return(NULL)
  }

  if (!inherits(taxon_evidence, "data.frame")) {
    stop("`taxon_evidence` must be a data.frame.", call. = FALSE)
  }

  required_columns <- .taxon_evidence_required_columns()
  missing_columns <- setdiff(required_columns, colnames(taxon_evidence))
  if (length(missing_columns) > 0) {
    stop(
      "`taxon_evidence` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  empty_columns <- required_columns[
    vapply(taxon_evidence[required_columns], .is_entirely_empty_column, logical(1))
  ]
  if (length(empty_columns) > 0) {
    stop(
      "Required `taxon_evidence` columns must not be entirely empty: ",
      paste(empty_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ordered_columns <- c(
    required_columns,
    intersect(.taxon_evidence_optional_columns(), colnames(taxon_evidence))
  )
  ordered_columns <- c(ordered_columns, setdiff(colnames(taxon_evidence), ordered_columns))

  taxon_evidence[, ordered_columns, drop = FALSE]
}

#' Extract a taxonomy table from supported inputs
#'
#' @param x A data frame, matrix, phyloseq taxonomy table, or phyloseq object.
#'
#' @return A data frame containing taxonomy.
#' @noRd
extract_tax_table <- function(x) {
  if (methods::is(x, "phyloseq")) {
    tax_table <- phyloseq::tax_table(x, errorIfNULL = FALSE)
    if (is.null(tax_table)) {
      stop("`x` does not contain a taxonomy table.", call. = FALSE)
    }
    tax_table <- as.data.frame(
      methods::as(tax_table, "matrix"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (inherits(x, "taxonomyTable")) {
    tax_table <- as.data.frame(
      methods::as(x, "matrix"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (inherits(x, "data.frame")) {
    tax_table <- as.data.frame(
      x,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (is.matrix(x)) {
    tax_table <- as.data.frame(
      x,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    stop(
      "`x` must be a data.frame taxonomy table or a phyloseq/speedyseq object.",
      call. = FALSE
    )
  }

  if (nrow(tax_table) < 1) {
    stop("`x` must contain at least one taxonomic row.", call. = FALSE)
  }
  if (ncol(tax_table) < 1) {
    stop("`x` must contain at least one taxonomy column.", call. = FALSE)
  }

  .add_taxon_id_column(tax_table)
}

#' Choose the most specific available taxon from a taxonomy row
#'
#' @param tax_row A named taxonomy row.
#' @param tax_ranks Character vector of ranks ordered from most specific to
#'   broadest.
#'
#' @return A list with `query_rank` and `query_name`.
#' @noRd
choose_query_taxon <- function(
  tax_row,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum", "kingdom")
) {
  tax_row <- .coerce_tax_row(tax_row)
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")

  for (rank in tax_ranks) {
    taxon_name <- .tax_rank_value(tax_row, rank)
    if (.is_informative_taxon(taxon_name)) {
      return(list(
        query_rank = tolower(rank),
        query_name = .clean_taxon_name(taxon_name)
      ))
    }
  }

  list(
    query_rank = NA_character_,
    query_name = NA_character_
  )
}

#' Build the flag_taxa prompt
#'
#' @param tax_table Data frame taxonomy table.
#' @inheritParams flag_taxa
#'
#' @return A character scalar prompt.
#' @noRd
build_flag_taxa_prompt <- function(
  tax_table,
  expected_environment,
  expected_habitat = NULL,
  expected_region = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum", "kingdom"),
  evidence_sources = c("worms", "obis", "gbif", "bold", "literature"),
  taxon_evidence = NULL
) {
  if (!inherits(tax_table, "data.frame")) {
    stop("`tax_table` must be a data.frame.", call. = FALSE)
  }

  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")
  evidence_sources <- .validate_character_vector(
    evidence_sources,
    "evidence_sources"
  )
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  prepared_taxonomy <- .prepare_flag_taxa_taxonomy(tax_table, tax_ranks)

  expected_habitat_text <- .format_optional_prompt_value(expected_habitat)
  expected_region_text <- .format_optional_prompt_value(expected_region)
  evidence_text <- paste(evidence_sources, collapse = ", ")
  json_fields <- .flag_taxa_json_output_fields("feature_id" %in% colnames(tax_table))
  user_evidence_section <- .format_taxon_evidence_prompt_section(taxon_evidence)

  paste(
    "You are reviewing molecular ecology/metabarcoding taxonomy for possible ecological or geographic incompatibilities.",
    "",
    "Task:",
    "Assess each row in the input taxonomy table against the expected environment, expected habitat, and expected region.",
    "Use the query_rank and query_name columns as the taxonomic group to assess, and use lineage for broader context.",
    "",
    "Expected context:",
    paste0("- expected_environment: ", expected_environment),
    paste0("- expected_habitat: ", expected_habitat_text),
    paste0("- expected_region: ", expected_region_text),
    paste0("- preferred evidence sources: ", evidence_text),
    user_evidence_section,
    "",
    "Return format:",
    "Return valid JSON only.",
    "Do not return markdown, a markdown table, a plain text table, code fences, or explanatory text outside the JSON.",
    "The JSON must be an array of objects, with one object for each input row in the same order.",
    "Each object must contain these keys:",
    paste(paste0("- ", json_fields), collapse = "\n"),
    "",
    "Allowed values:",
    "expected_environment_status:",
    "- compatible",
    "- incompatible",
    "- mixed_within_rank",
    "- unknown",
    "- insufficient_taxonomic_resolution",
    "",
    "expected_habitat_status:",
    "- compatible",
    "- incompatible",
    "- transient_or_allochthonous_possible",
    "- mixed_within_rank",
    "- unknown",
    "- insufficient_taxonomic_resolution",
    "",
    "expected_region_status:",
    "- known_in_region",
    "- known_near_region",
    "- known_elsewhere_only",
    "- restricted_elsewhere",
    "- no_distribution_evidence",
    "- unknown",
    "- not_assessed",
    "",
    "recommended_action:",
    "- retain",
    "- review",
    "- exclude",
    "",
    "Important biological rules:",
    "- Do not classify taxa as incompatible based on majority ecology.",
    "- A taxon is environmentally incompatible only when no known member of the queried taxonomic group is compatible with the expected environment.",
    "- If a broad taxonomic rank contains both compatible and incompatible organisms, use mixed_within_rank.",
    "- If taxonomy is too coarse to judge, use insufficient_taxonomic_resolution.",
    "- If evidence is missing, use unknown or no_distribution_evidence. Do not invent evidence.",
    "- Distinguish living habitat incompatibility from transient or allochthonous DNA.",
    "- Use exclude only for obvious incompatibilities.",
    "- Include references that justify categorisations.",
    "",
    "If expected_habitat is not specified, use unknown for expected_habitat_status unless evidence supports a more specific interpretation.",
    "If expected_region is not specified, use not_assessed for expected_region_status.",
    "Write concise rationales and include references for each row. Do not add free-text commentary outside the JSON.",
    "",
    "Input taxonomy table with query fields, tab-separated:",
    .format_prompt_table(prepared_taxonomy),
    sep = "\n"
  )
}

.prepare_flag_taxa_taxonomy <- function(tax_table, tax_ranks) {
  reserved_columns <- .flag_taxa_output_columns()
  collisions <- colnames(tax_table)[tolower(colnames(tax_table)) %in% reserved_columns]
  if (length(collisions) > 0) {
    stop(
      "`x` already contains columns reserved by `flag_taxa()`: ",
      paste(collisions, collapse = ", "),
      call. = FALSE
    )
  }

  matched_ranks <- tolower(tax_ranks) %in% tolower(colnames(tax_table))
  if (!any(matched_ranks)) {
    stop(
      "`x` must contain at least one requested taxonomic rank column. ",
      "Expected one of: ",
      paste(tax_ranks, collapse = ", "),
      call. = FALSE
    )
  }

  query_taxa <- lapply(seq_len(nrow(tax_table)), function(row_index) {
    choose_query_taxon(tax_table[row_index, , drop = FALSE], tax_ranks = tax_ranks)
  })

  tax_table$query_rank <- vapply(
    query_taxa,
    function(query_taxon) query_taxon$query_rank,
    character(1)
  )
  tax_table$query_name <- vapply(
    query_taxa,
    function(query_taxon) query_taxon$query_name,
    character(1)
  )
  tax_table$lineage <- vapply(
    seq_len(nrow(tax_table)),
    function(row_index) .build_lineage(tax_table[row_index, , drop = FALSE], tax_ranks),
    character(1)
  )

  tax_table
}

.flag_taxa_output_columns <- function() {
  c(
    "query_rank",
    "query_name",
    "lineage",
    "expected_environment",
    "expected_environment_status",
    "expected_habitat",
    "expected_habitat_status",
    "expected_region",
    "expected_region_status",
    "recommended_action",
    "rationale",
    "references"
  )
}

.flag_taxa_json_output_fields <- function(include_feature_id) {
  fields <- .flag_taxa_output_columns()
  if (isTRUE(include_feature_id)) {
    fields <- c("feature_id", fields)
  }

  fields
}

.taxon_evidence_required_columns <- function() {
  c(
    "taxon_name",
    "taxon_rank",
    "source",
    "evidence_type",
    "evidence_summary",
    "reference"
  )
}

.taxon_evidence_optional_columns <- function() {
  c(
    "accepted_name",
    "accepted_rank",
    "source_taxon_id",
    "source_record_id",
    "environment",
    "habitat",
    "region",
    "locality",
    "decimal_latitude",
    "decimal_longitude",
    "basis_of_record",
    "occurrence_count",
    "reference_url",
    "doi",
    "checked_at"
  )
}

.is_entirely_empty_column <- function(column) {
  values <- as.character(column)
  all(is.na(values) | !nzchar(trimws(values)))
}

.format_taxon_evidence_prompt_section <- function(taxon_evidence) {
  if (is.null(taxon_evidence)) {
    return("")
  }

  paste(
    "",
    "User-provided taxon evidence",
    "Use this evidence before performing online searches.",
    "Do not search online for taxa where the provided evidence is sufficient to assign the requested statuses.",
    "Use online sources only when the provided evidence is incomplete, missing, or conflicting.",
    "Evidence table, tab-separated:",
    .format_prompt_table(taxon_evidence),
    sep = "\n"
  )
}

.flag_taxa_allowed_values <- function() {
  list(
    expected_environment_status = c(
      "compatible",
      "incompatible",
      "mixed_within_rank",
      "unknown",
      "insufficient_taxonomic_resolution"
    ),
    expected_habitat_status = c(
      "compatible",
      "incompatible",
      "transient_or_allochthonous_possible",
      "mixed_within_rank",
      "unknown",
      "insufficient_taxonomic_resolution"
    ),
    expected_region_status = c(
      "known_in_region",
      "known_near_region",
      "known_elsewhere_only",
      "restricted_elsewhere",
      "no_distribution_evidence",
      "unknown",
      "not_assessed"
    ),
    recommended_action = c(
      "retain",
      "review",
      "exclude"
    )
  )
}

.validate_flag_taxa_enum_column <- function(result, column_name, allowed_values) {
  values <- as.character(result[[column_name]])
  invalid_values <- unique(values[is.na(values) | !(values %in% allowed_values)])

  if (length(invalid_values) > 0) {
    invalid_values[is.na(invalid_values)] <- "<NA>"
    stop(
      "`",
      column_name,
      "` contains invalid values: ",
      paste(invalid_values, collapse = ", "),
      ". Allowed values are: ",
      paste(allowed_values, collapse = ", "),
      call. = FALSE
    )
  }
}

.validate_flag_taxa_region_assessment <- function(result, expected_region) {
  region_status <- as.character(result$expected_region_status)
  not_assessed <- region_status == "not_assessed"

  if (is.null(expected_region) && any(!not_assessed)) {
    invalid_values <- unique(region_status[!not_assessed])
    stop(
      "`expected_region_status` must be `not_assessed` when ",
      "`expected_region` is NULL. Found: ",
      paste(invalid_values, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(expected_region) && any(not_assessed)) {
    stop(
      "`expected_region_status` must not be `not_assessed` when ",
      "`expected_region` is provided.",
      call. = FALSE
    )
  }
}

.add_taxon_id_column <- function(tax_table) {
  row_names <- rownames(tax_table)
  default_row_names <- as.character(seq_len(nrow(tax_table)))
  if (
    !is.null(row_names) &&
      !identical(row_names, default_row_names) &&
      !("taxon_id" %in% colnames(tax_table))
  ) {
    tax_table <- cbind(
      taxon_id = row_names,
      tax_table,
      stringsAsFactors = FALSE
    )
  }

  rownames(tax_table) <- NULL
  tax_table
}

.coerce_tax_row <- function(tax_row) {
  if (inherits(tax_row, "data.frame")) {
    if (nrow(tax_row) != 1) {
      stop("`tax_row` must contain exactly one row.", call. = FALSE)
    }
    out <- vapply(tax_row, function(value) as.character(value[[1]]), character(1))
  } else if (is.list(tax_row)) {
    out <- vapply(tax_row, function(value) as.character(value[[1]]), character(1))
  } else {
    out <- as.character(tax_row)
    names(out) <- names(tax_row)
  }

  if (is.null(names(out)) || any(names(out) == "")) {
    stop("`tax_row` must be a named taxonomy row.", call. = FALSE)
  }

  out
}

.tax_rank_value <- function(tax_row, rank) {
  rank_index <- which(tolower(names(tax_row)) == tolower(rank))
  if (length(rank_index) == 0) {
    return(NA_character_)
  }

  as.character(tax_row[[rank_index[[1]]]])
}

.build_lineage <- function(tax_row, tax_ranks) {
  tax_row <- .coerce_tax_row(tax_row)
  lineage_ranks <- rev(tax_ranks)
  lineage <- character()

  for (rank in lineage_ranks) {
    taxon_name <- .tax_rank_value(tax_row, rank)
    if (.is_informative_taxon(taxon_name)) {
      lineage <- c(
        lineage,
        paste0(tolower(rank), ": ", .clean_taxon_name(taxon_name))
      )
    }
  }

  if (length(lineage) == 0) {
    return(NA_character_)
  }

  paste(lineage, collapse = "; ")
}

.is_informative_taxon <- function(taxon_name) {
  if (length(taxon_name) == 0 || is.na(taxon_name[[1]])) {
    return(FALSE)
  }

  cleaned_name <- tolower(.clean_taxon_name(taxon_name[[1]]))
  unknown_values <- c(
    "",
    "na",
    "nan",
    "null",
    "none",
    "unknown",
    "unclassified",
    "unassigned",
    "undetermined",
    "unidentified",
    "not assigned"
  )

  !(cleaned_name %in% unknown_values)
}

.clean_taxon_name <- function(taxon_name) {
  taxon_name <- trimws(as.character(taxon_name[[1]]))
  taxon_name <- gsub("[\r\n\t]+", " ", taxon_name)
  taxon_name <- sub("^[[:alpha:]]__", "", taxon_name)
  trimws(taxon_name)
}

.format_prompt_table <- function(tax_table) {
  formatted <- as.data.frame(
    tax_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  formatted[] <- lapply(formatted, function(column) {
    column <- as.character(column)
    column[is.na(column)] <- ""
    gsub("[\r\n\t]+", " ", column)
  })

  rows <- apply(formatted, 1, paste, collapse = "\t")
  paste(c(paste(colnames(formatted), collapse = "\t"), rows), collapse = "\n")
}

.format_optional_prompt_value <- function(value) {
  if (is.null(value)) {
    return("not specified")
  }

  value
}

.validate_required_scalar_character <- function(value, name) {
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

.validate_optional_scalar_character <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }

  .validate_required_scalar_character(value, name)
}

.validate_character_vector <- function(value, name) {
  if (
    !is.character(value) ||
      length(value) < 1 ||
      any(is.na(value)) ||
      any(!nzchar(trimws(value)))
  ) {
    stop("`", name, "` must be a non-empty character vector.", call. = FALSE)
  }

  unique(tolower(trimws(value)))
}

.validate_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }

  value
}
