#' Deterministically assess taxa against local evidence
#'
#' @description
#' `assess_taxa_evidence()` runs the evidence-first, deterministic part of the
#' `flag_taxa()` workflow. It extracts the highest useful taxonomic assignment
#' for each input feature, matches user-supplied or locally read evidence, and
#' returns environment, habitat, region, ecological, occurrence, and recommended
#' action fields.
#'
#' Missing evidence, failed lookups, no WoRMS match, absent CCZ records, and
#' absent habitat or region records are treated as uncertainty rather than as
#' positive incompatibility evidence. OBIS checklist results are used as
#' regional occurrence evidence only: no OBIS records means missing regional
#' evidence, not evidence of absence, and OBIS is not treated as proof of
#' resident habitat use.
#'
#' `expected_habitat_status = "no_known_habitat_evidence"` means an expected
#' habitat was supplied, but deterministic evidence sources did not document the
#' taxon in that habitat. `expected_habitat_status = "unknown"` is reserved for
#' unassessed or unusable habitat context. Missing habitat evidence triggers
#' review rather than exclusion unless positive incompatibility evidence is
#' present. A retained row means available deterministic evidence supports
#' plausibility and no positive incompatibility was detected; it does not mean
#' every ecological detail is proven.
#'
#' Environment compatibility is assessed from trusted local evidence first:
#' user-provided taxon evidence, then WoRMS-derived environment flags when
#' requested. OBIS is primarily a regional occurrence layer. OBIS environment
#' flags are used only as fallback environment evidence when no matched
#' user/WoRMS environment evidence is available, and this fallback is noted in
#' the output. Terrestrial-only or freshwater-only trusted evidence in an
#' expected marine sample is treated as a clear environment incompatibility and
#' may be excluded deterministically. OBIS regional records do not override
#' trusted environment incompatibility.
#'
#' User-provided evidence and OBIS checklist evidence are kept as separate
#' evidence layers in the output. Both can contribute to
#' `expected_region_status` using the priority `known_in_region` >
#' `plausible_in_region` > `no_distribution_evidence`.
#'
#' FishBase/SeaLifeBase ecology evidence, when supplied through
#' `ecology_evidence` or `ecology_evidence_path`, is a habitat/ecology metadata
#' layer only. It may refine `expected_habitat_status`,
#' `ecological_status`, `occurrence_interpretation`,
#' `recommended_action`, `rationale`, and `ecology_evidence_*` output columns,
#' but it does not influence `expected_environment_status` or
#' `expected_region_status`.
#'
#' The function always returns a list. `res$summary` is a compact table for
#' reviewing taxa: taxon identity, lineage, key WoRMS/OBIS/ecology summaries,
#' then the deterministic status/action columns near the end. `res$details` is
#' the cleaned wide audit table. `res$worms`, `res$obis`, and `res$ecology`
#' contain source-specific evidence slices and are empty tibbles when no columns
#' are available for that source.
#'
#' The public ecology habitat field is `ecology_evidence_habitat_broad`, derived
#' from FishBase/SeaLifeBase broad habitat metadata. Redundant
#' `ecology_evidence_position` values are not shown in the user-facing outputs.
#' Aggregated semicolon-separated ecology terms are de-duplicated, and dominant
#' ecology fields in the detailed/audit tables contain a single value plus a
#' proportion. The summary keeps only compact ecology fields, including min/max
#' depth, and uses a short one-sentence rationale rather than the accumulated
#' audit rationale.
#' In `expected_habitat_type = "deep_sea_benthic"` context, Aves and Reptilia
#' are flagged as possible transient/allochthonous signals rather than resident
#' benthic organisms. Mammalia are not blanket-flagged by this rule because some
#' marine mammals dive deeply.
#'
#' @param x A taxonomy table as a data frame, matrix, phyloseq taxonomy table,
#'   or a phyloseq/speedyseq object containing a taxonomy table.
#' @param expected_environment Character scalar describing the expected broad
#'   environment, for example `"marine"`, `"freshwater"`, or `"terrestrial"`.
#' @param ... Reserved for backward compatibility with earlier calls.
#' @param expected_habitat_type Optional character scalar describing the broad
#'   habitat type used by deterministic ecology rules. Allowed values are
#'   `"unspecified"`, `"coastal_benthic"`, `"deep_sea_benthic"`,
#'   `"deep_sea_pelagic"`, `"pelagic"`, `"freshwater"`, `"terrestrial"`, and
#'   `"host_associated"`. `"deep_sea_benthic"` enables conservative treatment
#'   of birds and reptiles as likely transient/allochthonous in deep-sea
#'   benthic samples. If omitted or `NULL`, `"unspecified"` is used.
#' @param expected_region Optional character scalar describing the expected
#'   geographic region.
#' @param tax_ranks Character vector of taxonomic ranks to inspect, ordered
#'   from most specific to broadest.
#' @param evidence_sources Optional character vector of on-demand evidence
#'   sources to fetch or use before deterministic assessment. Currently
#'   `"worms"` and `"obis"` can trigger on-demand queries when their other
#'   inputs are supplied; `"ecology"` is accepted as an explicit local ecology
#'   evidence layer when `ecology_evidence` or `ecology_evidence_path` is
#'   supplied. Use `NULL` to avoid automatic evidence retrieval.
#' @param taxon_evidence Optional data frame containing user-supplied evidence
#'   about taxon environment, habitat, and/or region, including output from
#'   `fetch_worms_evidence()`.
#' @param taxon_evidence_path Optional path to a locally stored taxon evidence
#'   table. Standard `taxon_evidence` tables are read directly; CCZ-like WoRMS
#'   checklist files are standardised with `read_worms_ccz_evidence()`. This
#'   never downloads data.
#' @param ecology_evidence Optional FishBase/SeaLifeBase ecology evidence table,
#'   usually `build_fishbase_sealifebase_ecology_db()$combined`, or the full
#'   list returned by `build_fishbase_sealifebase_ecology_db()`. This evidence
#'   is used only for habitat/ecology interpretation and output metadata. It
#'   never changes `expected_environment_status` or
#'   `expected_region_status`.
#' @param ecology_evidence_path Optional path to a locally stored
#'   FishBase/SeaLifeBase ecology evidence table or RDS list containing a
#'   `combined` element. Supported formats are `.csv`, `.tsv`, `.txt`, `.rds`,
#'   `.xlsx`, and `.xls` when `readxl` is available.
#' @param obis_geometry Optional WKT geometry used for OBIS regional checklist
#'   queries. When `evidence_sources` includes `"obis"`, this function calls
#'   `robis::checklist(scientificname = query_name, geometry = obis_geometry)`
#'   for each unique query taxon. Raw OBIS occurrences, depth summaries, and
#'   caching are not used.
#' @param obis_region_label Optional label for `obis_geometry`. If `NULL`, the
#'   label defaults to `expected_region` when available, otherwise
#'   `"obis_geometry"`.
#' @param obis_region_relation Character scalar, either `"direct"` or
#'   `"supporting"`. Use `"direct"` when `obis_geometry` is the exact expected
#'   region and OBIS checklist records should support `known_in_region`. Use
#'   `"supporting"` when the geometry is broader or contextually relevant and
#'   records should support `plausible_in_region`.
#' @param taxon_evidence_region_relation Character scalar, either `"direct"` or
#'   `"supporting"`. Use `"direct"` when user-supplied `taxon_evidence` is from
#'   the exact expected region. Use `"supporting"` when user evidence is from a
#'   broader or related region.
#' @param taxon_evidence_region_label Optional label for the region or context
#'   represented by user-provided `taxon_evidence` or `taxon_evidence_path`.
#'   When `NULL`, the region is inferred from matched evidence rows when
#'   possible.
#' @param verbose Logical scalar. If `TRUE`, report progress messages.
#' @param worms_batch_size Positive integer scalar. Maximum number of unique
#'   taxa sent to WoRMS per live query batch when on-demand WoRMS evidence is
#'   requested. Defaults to `25`.
#' @param worms_max_tries Positive integer scalar. Maximum number of attempts
#'   for transient WoRMS/API failures per batch. Defaults to `3`.
#' @param worms_retry_sleep Non-negative numeric scalar. Seconds to wait between
#'   WoRMS retry attempts. Defaults to `5`.
#'
#' @return A list with `summary`, `details`, `worms`, `obis`, and `ecology`
#'   elements. `summary` is a compact tibble with one row per original input
#'   feature/taxon row and key columns such as `feature_id`, `taxon_name`,
#'   `taxon_rank`, `lineage`, `worms_environment`,
#'   `ecology_evidence_habitat_broad`, `ecology_evidence_depth_zone`,
#'   `ecology_evidence_depth_min`, `ecology_evidence_depth_max`,
#'   `obis_region_status`, deterministic status columns, `recommended_action`,
#'   short `rationale`, and `references`. It does not include internal
#'   `query_id`, redundant `ecology_evidence_position`, dominant ecology
#'   columns, or OBIS compact rank columns `obis_n_rank`/`obis_taxa_rank`.
#'   `details` is the cleaned wide audit table with the longer rationale and
#'   fuller evidence summaries. The source-specific tables contain the columns
#'   that are available for that run and are empty tibbles when no columns are
#'   available.
#' @export
#'
#' @examples
#' tax <- data.frame(
#'   phylum = "Chordata",
#'   family = "Salmonidae",
#'   genus = "Salmo",
#'   species = "Salmo salar"
#' )
#'
#' assess_taxa_evidence(
#'   tax,
#'   expected_environment = "marine",
#'   expected_habitat_type = "coastal_benthic"
#' )
#'
#' \dontrun{
#' # With a local FishBase/SeaLifeBase ecology DB and optional regional evidence:
#' res <- biohelper::assess_taxa_evidence(
#'   x = ps,
#'   expected_environment = "marine",
#'   expected_habitat_type = "deep_sea_benthic",
#'   expected_region = "Clarion-Clipperton Zone",
#'   evidence_sources = c("worms", "obis", "ecology"),
#'   ecology_evidence = fb_slb_db$combined
#' )
#'
#' res$summary
#' res$details
#' res$obis
#' res$ecology
#' }
assess_taxa_evidence <- function(
  x,
  expected_environment,
  ...,
  expected_habitat_type = c(
    "unspecified",
    "coastal_benthic",
    "deep_sea_benthic",
    "deep_sea_pelagic",
    "pelagic",
    "freshwater",
    "terrestrial",
    "host_associated"
  ),
  expected_region = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
  evidence_sources = NULL,
  taxon_evidence = NULL,
  taxon_evidence_path = NULL,
  ecology_evidence = NULL,
  ecology_evidence_path = NULL,
  obis_geometry = NULL,
  obis_region_label = NULL,
  obis_region_relation = c("direct", "supporting"),
  taxon_evidence_region_relation = c("direct", "supporting"),
  taxon_evidence_region_label = NULL,
  verbose = FALSE,
  worms_batch_size = 25,
  worms_max_tries = 3,
  worms_retry_sleep = 5
) {
  start_time <- Sys.time()
  worms_cache_path <- NULL
  worms_sleep <- 0.2
  continue_on_worms_error <- TRUE
  worms_verbose <- verbose
  max_evidence_rows_per_query <- .Machine$integer.max
  max_evidence_summary_chars <- 200L

  expected_environment <- .validate_required_scalar_character(
    expected_environment,
    "expected_environment"
  )
  dots <- .parse_assess_taxa_evidence_dots(...)
  expected_habitat <- dots$expected_habitat
  expected_habitat_type_arg <- if (missing(expected_habitat_type)) {
    NULL
  } else {
    expected_habitat_type
  }
  expected_habitat_type <- .validate_assess_expected_habitat_type(
    expected_habitat_type_arg,
    expected_habitat
  )
  if (is.null(expected_habitat)) {
    expected_habitat <- .assess_expected_habitat_from_type(expected_habitat_type)
  }
  expected_region <- .validate_optional_scalar_character(
    expected_region,
    "expected_region"
  )
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")
  evidence_sources <- .validate_assess_taxa_evidence_sources(evidence_sources)
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  taxon_evidence_path <- .validate_optional_path(taxon_evidence_path, "taxon_evidence_path")
  ecology_evidence <- .validate_assess_ecology_evidence(ecology_evidence)
  ecology_evidence_path <- .validate_optional_path(ecology_evidence_path, "ecology_evidence_path")
  obis_geometry <- .validate_optional_scalar_character(obis_geometry, "obis_geometry")
  obis_region_label <- .validate_optional_scalar_character(
    obis_region_label,
    "obis_region_label"
  )
  obis_region_relation <- match.arg(obis_region_relation)
  taxon_evidence_region_relation <- match.arg(taxon_evidence_region_relation)
  taxon_evidence_region_label <- .validate_optional_scalar_character(
    taxon_evidence_region_label,
    "taxon_evidence_region_label"
  )
  if (is.null(obis_region_label)) {
    obis_region_label <- if (!is.null(expected_region)) expected_region else "obis_geometry"
  }
  verbose <- .validate_logical_scalar(verbose, "verbose")
  worms_cache_path <- .validate_optional_cache_path(worms_cache_path)
  worms_sleep <- .validate_non_negative_numeric_scalar(worms_sleep, "worms_sleep")
  worms_batch_size <- .validate_positive_integer_scalar(worms_batch_size, "worms_batch_size")
  worms_max_tries <- .validate_positive_integer_scalar(worms_max_tries, "worms_max_tries")
  worms_retry_sleep <- .validate_non_negative_numeric_scalar(worms_retry_sleep, "worms_retry_sleep")
  worms_verbose <- .validate_logical_scalar(worms_verbose, "worms_verbose")
  max_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_evidence_summary_chars,
    "max_evidence_summary_chars"
  )

  tax_table <- extract_tax_table(x)
  .flag_taxa_progress(
    verbose,
    start_time,
    "input contains ",
    nrow(tax_table),
    " feature rows."
  )

  all_feature_map <- .prepare_flag_taxa_all_feature_map(
    tax_table = tax_table,
    tax_ranks = tax_ranks
  )
  feature_query_map <- .retained_flag_taxa_feature_query_map(all_feature_map)
  not_assessed_rows <- nrow(all_feature_map) - nrow(feature_query_map)
  .flag_taxa_progress(
    verbose,
    start_time,
    "retained ",
    nrow(feature_query_map),
    " feature rows with useful taxonomy; ",
    not_assessed_rows,
    " will be returned as not assessed."
  )

  unique_query_table <- .unique_flag_taxa_query_table(feature_query_map)
  feature_query_map <- .assign_flag_taxa_query_ids_to_features(
    feature_query_map,
    unique_query_table
  )
  .flag_taxa_progress(
    verbose,
    start_time,
    "assessing ",
    nrow(unique_query_table),
    " unique taxa from ",
    nrow(feature_query_map),
    " retained feature rows."
  )

  worms_evidence_sources <- intersect(evidence_sources, "worms")
  if (length(worms_evidence_sources) == 0) {
    worms_evidence_sources <- NULL
  }

  path_taxon_evidence <- .read_assess_taxon_evidence_path(taxon_evidence_path)
  user_taxon_evidence <- .combine_taxon_evidence_tables(taxon_evidence, path_taxon_evidence)
  path_ecology_evidence <- .read_assess_ecology_evidence_path(ecology_evidence_path)
  ecology_evidence <- .combine_assess_ecology_evidence_tables(ecology_evidence, path_ecology_evidence)
  taxon_evidence <- .collect_flag_taxa_requested_evidence(
    x = x,
    evidence_sources = worms_evidence_sources,
    taxon_evidence = user_taxon_evidence,
    ccz_evidence_path = NULL,
    worms_cache_path = worms_cache_path,
    worms_sleep = worms_sleep,
    worms_batch_size = worms_batch_size,
    worms_max_tries = worms_max_tries,
    worms_retry_sleep = worms_retry_sleep,
    continue_on_worms_error = continue_on_worms_error,
    worms_verbose = worms_verbose
  )

  obis_evidence <- NULL
  if ("obis" %in% evidence_sources) {
    if (is.null(obis_geometry)) {
      warning(
        "OBIS regional evidence requires `obis_geometry`; skipping OBIS checklist queries.",
        call. = FALSE
      )
    } else {
      obis_evidence <- fetch_obis_checklist_evidence(
        query_taxa = unique_query_table,
        obis_geometry = obis_geometry,
        obis_region_label = obis_region_label,
        obis_region_relation = obis_region_relation,
        verbose = verbose
      )
    }
  }

  evidence_stats <- .flag_taxa_evidence_match_stats(
    unique_query_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  .flag_taxa_progress(
    verbose,
    start_time,
    "matched ",
    evidence_stats$evidence_rows,
    " local evidence rows covering ",
    evidence_stats$query_taxa_with_evidence,
    " / ",
    nrow(unique_query_table),
    " unique taxa."
  )

  preliminary_judgement <- .build_flag_taxa_preliminary_judgement(
    unique_query_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    tax_ranks = tax_ranks
  )

  out <- .build_flag_taxa_evidence_first_output(
    preliminary_judgement = preliminary_judgement,
    feature_query_map = feature_query_map,
    unique_query_table = unique_query_table,
    selected_query_table = unique_query_table[0, , drop = FALSE],
    all_feature_map = all_feature_map,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    taxon_evidence = taxon_evidence,
    tool_evidence = NULL,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  out <- .apply_assess_taxa_regional_evidence(
    result = out,
    unique_query_table = unique_query_table,
    user_taxon_evidence = user_taxon_evidence,
    obis_evidence = obis_evidence,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    taxon_evidence_region_relation = taxon_evidence_region_relation,
    taxon_evidence_region_label = taxon_evidence_region_label,
    obis_region_relation = obis_region_relation
  )
  out <- .apply_assess_ecology_evidence(
    result = out,
    unique_query_table = unique_query_table,
    ecology_evidence = ecology_evidence,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_habitat_type = expected_habitat_type,
    expected_region = expected_region
  )
  out <- .apply_assess_deep_sea_benthic_vertebrate_rule(
    result = out,
    expected_habitat_type = expected_habitat_type
  )
  details <- .format_assess_taxa_evidence_output(out)
  summary <- .format_assess_taxa_evidence_summary(details)

  .flag_taxa_progress(
    verbose,
    start_time,
    "completed deterministic evidence assessment. Returned ",
    nrow(out),
    " feature-level rows.",
    total = TRUE
  )

  .format_assess_taxa_evidence_list(summary, details)
}

.format_assess_taxa_evidence_output <- function(result) {
  if (!inherits(result, "data.frame") || nrow(result) == 0) {
    return(result)
  }

  audit_columns <- intersect(
    c(
      "query_name",
      "query_rank",
      "matched_name",
      "matched_rank",
      "llm_prompt_selected",
      "llm_selected",
      "llm_tool_used",
      "llm_tool_name",
      "llm_tool_query",
      "llm_tool_evidence_summary",
      "llm_tool_references",
      "query_id",
      "ecology_evidence_position"
    ),
    colnames(result)
  )
  result <- result[, setdiff(colnames(result), audit_columns), drop = FALSE]
  result <- .add_assess_compact_rank_columns(result)

  leading <- intersect(
    c(
      "feature_id",
      "taxon_rank",
      "taxon_name",
      "lineage",
      "worms_environment",
      "expected_environment_status",
      "expected_habitat_status",
      "expected_region_status",
      "ecological_status",
      "occurrence_interpretation",
      "recommended_action",
      "rationale",
      "references",
      "evidence_basis",
      "evidence_sources",
      "evidence_summary",
      .assess_taxa_user_evidence_output_columns(),
      .assess_taxa_obis_output_columns(),
      .assess_taxa_ecology_output_columns(),
      "expected_environment",
      "expected_habitat",
      "expected_region"
    ),
    colnames(result)
  )
  trailing <- setdiff(colnames(result), leading)
  result[, c(leading, trailing), drop = FALSE]
}

.format_assess_taxa_evidence_summary <- function(details) {
  if (!inherits(details, "data.frame")) {
    return(details)
  }
  details <- .add_assess_compact_rank_columns(details)
  if ("rationale" %in% colnames(details)) {
    details$rationale <- .format_assess_short_rationale(details)
  }
  preferred <- .assess_taxa_summary_columns()
  out <- .assess_select_with_missing(details, preferred)
  out
}

.format_assess_taxa_evidence_list <- function(summary, details) {
  list(
    summary = summary,
    details = details,
    worms = .format_assess_taxa_worms_table(details),
    obis = .format_assess_taxa_obis_table(details),
    ecology = .format_assess_taxa_ecology_table(details)
  )
}

.assess_taxa_summary_columns <- function() {
  c(
    "feature_id",
    "taxon_name",
    "ecology_evidence_common_name",
    "taxon_rank",
    "lineage",
    "worms_environment",
    "ecology_evidence_habitat_broad",
    "ecology_evidence_depth_zone",
    "ecology_evidence_depth_min",
    "ecology_evidence_depth_max",
    "obis_region_status",
    "obis_taxa_returned",
    "obis_total_records_including_descendants",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale",
    "references"
  )
}

.format_assess_short_rationale <- function(details) {
  if (!inherits(details, "data.frame") || nrow(details) == 0) {
    return(character())
  }
  classes <- .assess_deep_sea_benthic_air_breathing_class(details)
  vapply(seq_len(nrow(details)), function(row_index) {
    env_status <- .assess_row_value(details, "expected_environment_status", row_index)
    habitat_status <- .assess_row_value(details, "expected_habitat_status", row_index)
    region_status <- .assess_row_value(details, "expected_region_status", row_index)
    action <- .assess_row_value(details, "recommended_action", row_index)
    occurrence <- .assess_row_value(details, "occurrence_interpretation", row_index)
    habitat <- .assess_ecology_split_values(.assess_row_value(details, "ecology_evidence_habitat_broad", row_index))
    evidence_sources <- .assess_row_value(details, "evidence_sources", row_index)

    if (identical(env_status, "incompatible")) {
      return("Freshwater/terrestrial evidence conflicts with expected sample.")
    }
    if (classes[[row_index]] %in% c("aves", "reptilia") &&
        identical(action, "flag_possible_exclusion") &&
        identical(occurrence, "possible_transient_or_allochthonous")) {
      label <- if (identical(classes[[row_index]], "aves")) "seabird" else "marine reptile"
      return(paste0("Marine ", label, " in deep-sea benthic sample; likely allochthonous."))
    }
    if (identical(occurrence, "possible_transient_or_allochthonous") &&
        any(habitat %in% "pelagic")) {
      return("Marine pelagic taxon in benthic sample; likely transient/allochthonous.")
    }
    if (identical(habitat_status, "compatible") &&
        any(habitat %in% c("benthic", "demersal", "benthopelagic"))) {
      return("Benthic/demersal ecology supports expected habitat.")
    }
    if (env_status %in% c("compatible", "mixed_within_rank", "possible_likely") &&
        region_status %in% c("known_in_region", "known_near_region", "plausible_in_region") &&
        habitat_status %in% c("unknown", "no_known_habitat_evidence")) {
      return("Marine and regionally plausible; habitat unknown.")
    }
    if (identical(action, "retain") &&
        (!.is_non_empty_value(evidence_sources) || identical(evidence_sources, "none"))) {
      return("No usable ecology evidence; retained by environment/region evidence.")
    }
    paste0(
      "Environment ", .format_summary_value(env_status),
      "; habitat ", .format_summary_value(habitat_status),
      "; region ", .format_summary_value(region_status),
      "."
    )
  }, character(1))
}

.assess_row_value <- function(x, column, row_index) {
  if (!(column %in% colnames(x))) {
    return(NA_character_)
  }
  as.character(x[[column]][[row_index]])
}

.assess_select_with_missing <- function(x, columns) {
  for (column in columns) {
    if (!(column %in% colnames(x))) {
      x[[column]] <- NA
    }
  }
  tibble::as_tibble(x[, columns, drop = FALSE])
}

.add_assess_compact_rank_columns <- function(result) {
  if (!inherits(result, "data.frame") || nrow(result) == 0) {
    return(result)
  }
  if (!all(c("obis_n_rank", "obis_taxa_rank") %in% colnames(result))) {
    obis <- .collapse_assess_obis_rank_columns(result)
    result$obis_n_rank <- obis$n
    result$obis_taxa_rank <- obis$taxa
  }
  if (!all(c("taxon_evidence_n_rank", "taxon_evidence_taxa_rank") %in% colnames(result))) {
    taxon <- .collapse_assess_taxon_evidence_rank_columns(result)
    result$taxon_evidence_n_rank <- taxon$n
    result$taxon_evidence_taxa_rank <- taxon$taxa
  }
  result
}

.collapse_assess_obis_rank_columns <- function(result) {
  .collapse_assess_rank_columns(
    result = result,
    n_prefix = "obis_n_",
    taxa_prefix = "obis_taxa_"
  )
}

.collapse_assess_taxon_evidence_rank_columns <- function(result) {
  .collapse_assess_rank_columns(
    result = result,
    n_prefix = "taxon_evidence_n_",
    taxa_prefix = "taxon_evidence_taxa_"
  )
}

.collapse_assess_rank_columns <- function(result, n_prefix, taxa_prefix) {
  ranks <- .assess_taxa_evidence_ranks()
  n <- rep(NA_integer_, nrow(result))
  taxa <- rep(NA_character_, nrow(result))
  for (row_index in seq_len(nrow(result))) {
    taxon_rank <- if ("taxon_rank" %in% colnames(result)) {
      tolower(as.character(result$taxon_rank[[row_index]]))
    } else {
      NA_character_
    }
    choice <- .assess_best_rank_column_for_row(result[row_index, , drop = FALSE], taxon_rank, n_prefix, taxa_prefix, ranks)
    n[[row_index]] <- choice$n
    taxa[[row_index]] <- choice$taxa
  }
  list(n = n, taxa = taxa)
}

.assess_best_rank_column_for_row <- function(row, taxon_rank, n_prefix, taxa_prefix, ranks) {
  rank_index <- match(taxon_rank, ranks)
  if (!is.na(rank_index)) {
    descendant_ranks <- if (rank_index < length(ranks)) {
      ranks[seq.int(rank_index + 1L, length(ranks))]
    } else {
      character()
    }
    candidate_ranks <- c(ranks[rank_index], descendant_ranks)
  } else {
    candidate_ranks <- rev(ranks)
  }
  candidate_ranks <- unique(candidate_ranks[!is.na(candidate_ranks)])
  fallback_ranks <- rev(ranks)
  for (rank in c(candidate_ranks, fallback_ranks)) {
    n_col <- paste0(n_prefix, rank)
    taxa_col <- paste0(taxa_prefix, rank)
    if (!(n_col %in% colnames(row))) {
      next
    }
    n <- suppressWarnings(as.integer(row[[n_col]][[1]]))
    taxa <- if (taxa_col %in% colnames(row)) as.character(row[[taxa_col]][[1]]) else NA_character_
    if (!is.na(n) && n > 0) {
      return(list(n = n, taxa = taxa))
    }
  }
  list(n = NA_integer_, taxa = NA_character_)
}

.format_assess_taxa_worms_table <- function(details) {
  .assess_select_existing(
    details,
    c("feature_id", "taxon_name", "taxon_rank", "worms_environment", "expected_environment_status", "references")
  )
}

.format_assess_taxa_obis_table <- function(details) {
  .assess_select_existing(
    .add_assess_compact_rank_columns(details),
    c(
      "feature_id",
      "taxon_name",
      "taxon_rank",
      .assess_taxa_obis_output_columns(),
      "obis_n_rank",
      "obis_taxa_rank"
    )
  )
}

.format_assess_taxa_ecology_table <- function(details) {
  .assess_select_existing(
    details,
    c(
      "feature_id",
      "taxon_name",
      "taxon_rank",
      .assess_taxa_ecology_output_columns()
    )
  )
}

.assess_select_existing <- function(x, columns) {
  columns <- intersect(columns, colnames(x))
  if (length(columns) == 0) {
    return(tibble::as_tibble(data.frame()))
  }
  tibble::as_tibble(x[, columns, drop = FALSE])
}

.validate_assess_taxa_evidence_sources <- function(evidence_sources) {
  if (is.null(evidence_sources)) {
    return(NULL)
  }

  evidence_sources <- .validate_character_vector(
    evidence_sources,
    "evidence_sources"
  )
  supported_sources <- c("worms", "obis", "ecology")
  unsupported_sources <- setdiff(evidence_sources, supported_sources)
  if (length(unsupported_sources) > 0) {
    stop(
      "`evidence_sources` contains unsupported value(s): ",
      paste(unsupported_sources, collapse = ", "),
      ". Currently supported values are: worms, obis, ecology.",
      call. = FALSE
    )
  }

  evidence_sources
}

.parse_assess_taxa_evidence_dots <- function(...) {
  dots <- list(...)
  if (length(dots) == 0) {
    return(list(expected_habitat = NULL))
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }
  unnamed <- !nzchar(dot_names)
  if (any(unnamed)) {
    stop(
      "Unused positional argument(s). Use `expected_habitat_type` for broad habitat context.",
      call. = FALSE
    )
  }

  unexpected <- setdiff(dot_names, "expected_habitat")
  if (length(unexpected) > 0) {
    stop(
      "Unused argument(s): ",
      paste(unexpected, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (sum(dot_names == "expected_habitat") > 1L) {
    stop("`expected_habitat` was supplied more than once.", call. = FALSE)
  }

  list(
    expected_habitat = .validate_optional_scalar_character(
      dots$expected_habitat,
      "expected_habitat"
    )
  )
}

.validate_assess_expected_habitat_type <- function(expected_habitat_type, expected_habitat) {
  if (is.null(expected_habitat_type)) {
    return(.infer_assess_expected_habitat_type(expected_habitat))
  }
  expected_habitat_type <- .validate_required_scalar_character(
    expected_habitat_type,
    "expected_habitat_type"
  )
  expected_habitat_type <- tolower(trimws(expected_habitat_type))
  allowed <- .assess_expected_habitat_types()
  if (!(expected_habitat_type %in% allowed)) {
    stop(
      "`expected_habitat_type` must be one of: ",
      paste(allowed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  expected_habitat_type
}

.assess_expected_habitat_from_type <- function(expected_habitat_type) {
  expected_habitat_type <- tolower(trimws(as.character(expected_habitat_type[[1]])))
  switch(
    expected_habitat_type,
    coastal_benthic = "coastal",
    deep_sea_benthic = "deep sea benthic",
    deep_sea_pelagic = "deep sea pelagic",
    pelagic = "pelagic",
    freshwater = "freshwater",
    terrestrial = "terrestrial",
    host_associated = "host associated",
    unspecified = NULL,
    NULL
  )
}

.assess_expected_habitat_types <- function() {
  c(
    "unspecified",
    "coastal_benthic",
    "deep_sea_benthic",
    "deep_sea_pelagic",
    "pelagic",
    "freshwater",
    "terrestrial",
    "host_associated"
  )
}

.infer_assess_expected_habitat_type <- function(expected_habitat) {
  if (is.null(expected_habitat) || !.is_non_empty_value(expected_habitat)) {
    return("unspecified")
  }
  habitat <- tolower(as.character(expected_habitat))
  if (grepl("freshwater|fresh water|river|stream|lake|pond", habitat)) {
    return("freshwater")
  }
  if (grepl("terrestrial|soil|land|forest|grassland|desert", habitat)) {
    return("terrestrial")
  }
  if (grepl("host|parasite|symbiont|epiphyt", habitat)) {
    return("host_associated")
  }
  is_deep <- grepl("deep[ -]?sea|abyss|hadal|bathyal|aphotic", habitat)
  is_benthic <- grepl(
    "benthic|benthos|sediment|seafloor|bottom|demersal|nodule|polymetallic",
    habitat
  )
  is_pelagic <- grepl("pelagic|water column|plankton|nekton", habitat)
  is_coastal <- grepl("coastal|shore|intertidal|subtidal|reef|kelp|estuar", habitat)
  if (isTRUE(is_deep && is_pelagic && !is_benthic)) {
    return("deep_sea_pelagic")
  }
  if (isTRUE(is_deep)) {
    return("deep_sea_benthic")
  }
  if (isTRUE(is_coastal && (is_benthic || grepl("reef|kelp|intertidal|subtidal|estuar", habitat)))) {
    return("coastal_benthic")
  }
  if (isTRUE(is_pelagic)) {
    return("pelagic")
  }
  "unspecified"
}

.read_assess_taxon_evidence_path <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }
  raw <- .read_assess_taxon_evidence_file(path)
  standard <- tryCatch(
    validate_taxon_evidence(raw),
    error = function(error) NULL
  )
  if (!is.null(standard)) {
    return(standard)
  }
  .flag_taxa_read_worms_ccz_evidence(path)
}

.read_assess_taxon_evidence_file <- function(path) {
  path <- path.expand(path)
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "rds")) {
    return(readRDS(path))
  }
  if (extension %in% c("csv")) {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (extension %in% c("tsv", "txt")) {
    return(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (extension %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Reading Excel taxon evidence files requires the `readxl` package. ",
        "Install it with install.packages(\"readxl\") or provide a CSV/TSV/RDS file.",
        call. = FALSE
      )
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  stop(
    "`taxon_evidence_path` must point to a .csv, .tsv, .txt, .rds, .xlsx, or .xls file.",
    call. = FALSE
  )
}

.read_assess_ecology_evidence_path <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }
  .validate_assess_ecology_evidence(.read_assess_taxon_evidence_file(path))
}

.validate_assess_ecology_evidence <- function(ecology_evidence) {
  if (is.null(ecology_evidence)) {
    return(NULL)
  }
  if (is.list(ecology_evidence) &&
      !inherits(ecology_evidence, "data.frame") &&
      "combined" %in% names(ecology_evidence)) {
    ecology_evidence <- ecology_evidence$combined
  }
  if (!inherits(ecology_evidence, "data.frame")) {
    stop(
      "`ecology_evidence` must be a data.frame or the list returned by ",
      "`build_fishbase_sealifebase_ecology_db()`.",
      call. = FALSE
    )
  }

  out <- janitor::clean_names(as.data.frame(ecology_evidence, stringsAsFactors = FALSE))
  required <- c("taxon_name", "taxon_rank")
  missing <- setdiff(required, colnames(out))
  if (length(missing) > 0) {
    stop(
      "`ecology_evidence` is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!("source" %in% colnames(out))) {
    out$source <- "fishbase_sealifebase"
  }

  out$taxon_rank <- tolower(trimws(as.character(out$taxon_rank)))
  standard_ranks <- .assess_taxa_ecology_standard_ranks()
  out <- out[out$taxon_rank %in% standard_ranks, , drop = FALSE]

  for (column_name in .assess_taxa_ecology_known_columns()) {
    if (!(column_name %in% colnames(out))) {
      out[[column_name]] <- NA
    }
  }

  tibble::as_tibble(out)
}

.combine_assess_ecology_evidence_tables <- function(...) {
  tables <- Filter(Negate(is.null), list(...))
  if (length(tables) == 0) {
    return(NULL)
  }
  .validate_assess_ecology_evidence(dplyr::bind_rows(tables))
}

.assess_taxa_ecology_standard_ranks <- function() {
  c("kingdom", "phylum", "class", "order", "family", "genus", "species")
}

.assess_taxa_ecology_known_columns <- function() {
  c(
    "source",
    "spec_code",
    "taxon_rank",
    "taxon_name",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "common_name",
    "rank_water_salinity",
    "rank_body_shape",
    "rank_habitat_note",
    "rank_distribution",
    "rank_diagnosis",
    "rank_remarks",
    "rank_etymology",
    "ecology_environment",
    "ecology_position",
    "ecology_habitat_broad",
    "ecology_depth_zone",
    "ecology_zone",
    "ecology_water_column_zone",
    "ecology_mobility",
    "ecology_size_class",
    "ecology_substrate",
    "ecology_special_habitat",
    "ecology_comments",
    "n_species",
    "n_species_with_habitat_broad",
    "dominant_habitat_broad",
    "dominant_habitat_broad_prop",
    "dominant_depth_zone",
    "dominant_depth_zone_prop",
    "depth_min",
    "depth_max",
    "common_depth_min",
    "common_depth_max",
    "ecology_evidence_strength"
  )
}

summarise_assess_taxon_evidence_matches <- function(
  unique_query_table,
  taxon_evidence,
  taxon_evidence_region_label = NULL,
  taxon_evidence_region_relation = c("direct", "supporting")
) {
  taxon_evidence_region_relation <- match.arg(taxon_evidence_region_relation)
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  if (is.null(taxon_evidence)) {
    return(NULL)
  }
  if (!inherits(unique_query_table, "data.frame") ||
      !all(c("query_id", "query_name", "query_rank") %in% colnames(unique_query_table))) {
    stop(
      "`unique_query_table` must contain `query_id`, `query_name`, and `query_rank`.",
      call. = FALSE
    )
  }
  if (nrow(unique_query_table) == 0) {
    return(.empty_assess_taxon_evidence_summary())
  }

  rows <- lapply(seq_len(nrow(unique_query_table)), function(row_index) {
    evidence_rows <- .raw_flag_taxa_relevant_evidence(
      query_table = unique_query_table,
      taxon_evidence = taxon_evidence,
      row_index = row_index
    )
    .summarise_assess_taxon_evidence_rows(
      evidence_rows = evidence_rows,
      query_id = as.character(unique_query_table$query_id[[row_index]]),
      region_label = taxon_evidence_region_label,
      region_relation = taxon_evidence_region_relation
    )
  })

  out <- .bind_rows_fill(rows)
  rownames(out) <- NULL
  out
}

.summarise_assess_taxon_evidence_rows <- function(
  evidence_rows,
  query_id,
  region_label,
  region_relation
) {
  out <- .empty_assess_taxon_evidence_summary(query_id = query_id)
  out$taxon_evidence_region_relation <- region_relation

  if (is.null(evidence_rows) || nrow(evidence_rows) == 0) {
    return(out)
  }

  evidence_names <- as.character(evidence_rows$taxon_name)
  evidence_names <- unique(evidence_names[.is_non_empty_value(evidence_names)])
  out$taxon_evidence_taxa_returned <- length(evidence_names)
  out$taxon_evidence_environment <- .collapse_distinct_values(
    .flag_taxa_evidence_values(evidence_rows, "environment"),
    max_values = 8
  )
  out$taxon_evidence_region <- .infer_assess_taxon_evidence_region(
    evidence_rows = evidence_rows,
    region_label = region_label
  )

  rank_values <- tolower(as.character(evidence_rows$taxon_rank))
  taxon_names <- as.character(evidence_rows$taxon_name)
  for (rank in .assess_taxa_evidence_ranks()) {
    rank_summary <- .taxon_evidence_rank_summary(
      rank = rank,
      rank_values = rank_values,
      taxon_names = taxon_names
    )
    out[[paste0("taxon_evidence_n_", rank)]] <- rank_summary$n
    out[[paste0("taxon_evidence_taxa_", rank)]] <- rank_summary$taxa
  }

  out
}

.empty_assess_taxon_evidence_summary <- function(query_id = character()) {
  out <- data.frame(
    query_id = as.character(query_id),
    taxon_evidence_region = NA_character_,
    taxon_evidence_region_relation = NA_character_,
    taxon_evidence_taxa_returned = 0L,
    taxon_evidence_environment = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (rank in .assess_taxa_evidence_ranks()) {
    out[[paste0("taxon_evidence_n_", rank)]] <- 0L
    out[[paste0("taxon_evidence_taxa_", rank)]] <- NA_character_
  }
  out
}

.assess_taxa_evidence_ranks <- function() {
  c("phylum", "class", "order", "family", "genus", "species")
}

.assess_taxa_user_evidence_output_columns <- function() {
  unlist(
    c(
      "taxon_evidence_region",
      "taxon_evidence_region_relation",
      "taxon_evidence_taxa_returned",
      "taxon_evidence_environment",
      lapply(.assess_taxa_evidence_ranks(), function(rank) {
        c(paste0("taxon_evidence_n_", rank), paste0("taxon_evidence_taxa_", rank))
      })
    ),
    use.names = FALSE
  )
}

summarise_assess_ecology_evidence_matches <- function(unique_query_table, ecology_evidence) {
  ecology_evidence <- .validate_assess_ecology_evidence(ecology_evidence)
  if (is.null(ecology_evidence)) {
    return(NULL)
  }
  if (!inherits(unique_query_table, "data.frame") ||
      !all(c("query_id", "query_name", "query_rank") %in% colnames(unique_query_table))) {
    stop(
      "`unique_query_table` must contain `query_id`, `query_name`, and `query_rank`.",
      call. = FALSE
    )
  }
  if (nrow(unique_query_table) == 0) {
    return(.empty_assess_ecology_evidence_summary())
  }

  rows <- lapply(seq_len(nrow(unique_query_table)), function(row_index) {
    matched <- .assess_ecology_match_indices(
      query_table = unique_query_table,
      ecology_evidence = ecology_evidence,
      row_index = row_index
    )
    evidence_rows <- ecology_evidence[matched, , drop = FALSE]
    .summarise_assess_ecology_evidence_rows(
      evidence_rows = evidence_rows,
      query_id = as.character(unique_query_table$query_id[[row_index]]),
      query_name = as.character(unique_query_table$query_name[[row_index]]),
      query_rank = as.character(unique_query_table$query_rank[[row_index]])
    )
  })

  out <- .bind_rows_fill(rows)
  rownames(out) <- NULL
  out
}

.assess_ecology_match_indices <- function(query_table, ecology_evidence, row_index) {
  query_name <- query_table$query_name[[row_index]]
  query_rank <- tolower(as.character(query_table$query_rank[[row_index]]))
  query_key <- .taxon_match_key(query_name)
  matched <- integer()

  name_keys <- .taxon_match_key(ecology_evidence$taxon_name)
  matched <- c(matched, which(name_keys == query_key))

  lineage_columns <- c(
    kingdom = "kingdom",
    phylum = "phylum",
    class = "class",
    order = "order",
    family = "family",
    genus = "genus"
  )
  if (query_rank %in% names(lineage_columns)) {
    lineage_column <- lineage_columns[[query_rank]]
    if (lineage_column %in% colnames(ecology_evidence)) {
      lineage_keys <- .taxon_match_key(ecology_evidence[[lineage_column]])
      matched <- c(matched, which(lineage_keys == query_key))
    }
  }

  unique(matched)
}

.summarise_assess_ecology_evidence_rows <- function(
  evidence_rows,
  query_id,
  query_name = NA_character_,
  query_rank = NA_character_
) {
  out <- .empty_assess_ecology_evidence_summary(query_id = query_id)
  if (is.null(evidence_rows) || nrow(evidence_rows) == 0) {
    return(out)
  }

  out$ecology_evidence_sources <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "source"),
    max_values = 6
  )
  out$ecology_evidence_taxa_returned <- length(unique(
    as.character(evidence_rows$taxon_name)[.is_non_empty_value(evidence_rows$taxon_name)]
  ))
  out$ecology_evidence_taxa <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "taxon_name"),
    max_values = 8
  )
  out$ecology_evidence_rank <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "taxon_rank"),
    max_values = 7
  )
  out$ecology_evidence_environment <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "ecology_environment"),
    max_values = 8
  )
  out$ecology_evidence_position <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "ecology_position"),
    max_values = 8
  )
  out$ecology_evidence_habitat_broad <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "ecology_habitat_broad"),
    max_values = 8
  )
  out$ecology_evidence_depth_zone <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "ecology_depth_zone"),
    max_values = 8
  )
  habitat_dom <- .assess_dominant_semicolon_value(
    .assess_ecology_column(evidence_rows, "ecology_habitat_broad"),
    fallback = .assess_ecology_column(evidence_rows, "dominant_habitat_broad")
  )
  out$ecology_evidence_dominant_habitat_broad <- habitat_dom$value
  out$ecology_evidence_dominant_habitat_broad_prop <- habitat_dom$prop
  depth_dom <- .assess_dominant_semicolon_value(
    .assess_ecology_column(evidence_rows, "ecology_depth_zone"),
    fallback = .assess_ecology_column(evidence_rows, "dominant_depth_zone")
  )
  out$ecology_evidence_dominant_depth_zone <- depth_dom$value
  out$ecology_evidence_dominant_depth_zone_prop <- depth_dom$prop
  out$ecology_evidence_strength <- .collapse_distinct_values(
    .assess_ecology_column(evidence_rows, "ecology_evidence_strength"),
    max_values = 5
  )
  out$ecology_evidence_n_species <- .assess_ecology_sum_numeric(
    .assess_ecology_column(evidence_rows, "n_species")
  )
  out$ecology_evidence_common_name <- .select_assess_ecology_common_name(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank
  )
  out$ecology_evidence_rank_habitat_note <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_habitat_note"
  )
  out$ecology_evidence_rank_water_salinity <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_water_salinity"
  )
  out$ecology_evidence_rank_body_shape <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_body_shape"
  )
  out$ecology_evidence_rank_distribution <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_distribution"
  )
  out$ecology_evidence_rank_diagnosis <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_diagnosis"
  )
  out$ecology_evidence_rank_remarks <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_remarks"
  )
  out$ecology_evidence_rank_etymology <- .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "rank_etymology"
  )
  out$ecology_evidence_depth_min <- .assess_ecology_min_numeric(
    .assess_ecology_column(evidence_rows, "depth_min")
  )
  out$ecology_evidence_depth_max <- .assess_ecology_max_numeric(
    .assess_ecology_column(evidence_rows, "depth_max")
  )
  out$ecology_evidence_common_depth_min <- .assess_ecology_min_numeric(
    .assess_ecology_column(evidence_rows, "common_depth_min")
  )
  out$ecology_evidence_common_depth_max <- .assess_ecology_max_numeric(
    .assess_ecology_column(evidence_rows, "common_depth_max")
  )
  out$ecology_evidence_summary <- .format_assess_ecology_summary(out)
  out
}

.empty_assess_ecology_evidence_summary <- function(query_id = character()) {
  data.frame(
    query_id = as.character(query_id),
    ecology_evidence_sources = NA_character_,
    ecology_evidence_taxa_returned = 0L,
    ecology_evidence_taxa = NA_character_,
    ecology_evidence_rank = NA_character_,
    ecology_evidence_environment = NA_character_,
    ecology_evidence_position = NA_character_,
    ecology_evidence_habitat_broad = NA_character_,
    ecology_evidence_depth_zone = NA_character_,
    ecology_evidence_dominant_habitat_broad = NA_character_,
    ecology_evidence_dominant_habitat_broad_prop = NA_real_,
    ecology_evidence_dominant_depth_zone = NA_character_,
    ecology_evidence_dominant_depth_zone_prop = NA_real_,
    ecology_evidence_strength = NA_character_,
    ecology_evidence_n_species = NA_real_,
    ecology_evidence_common_name = NA_character_,
    ecology_evidence_rank_habitat_note = NA_character_,
    ecology_evidence_rank_water_salinity = NA_character_,
    ecology_evidence_rank_body_shape = NA_character_,
    ecology_evidence_rank_distribution = NA_character_,
    ecology_evidence_rank_diagnosis = NA_character_,
    ecology_evidence_rank_remarks = NA_character_,
    ecology_evidence_rank_etymology = NA_character_,
    ecology_evidence_depth_min = NA_real_,
    ecology_evidence_depth_max = NA_real_,
    ecology_evidence_common_depth_min = NA_real_,
    ecology_evidence_common_depth_max = NA_real_,
    ecology_evidence_summary = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.assess_taxa_ecology_output_columns <- function() {
  colnames(.empty_assess_ecology_evidence_summary(query_id = NA_character_))[
    colnames(.empty_assess_ecology_evidence_summary(query_id = NA_character_)) != "query_id"
  ]
}

.assess_ecology_column <- function(evidence_rows, column_name) {
  if (!(column_name %in% colnames(evidence_rows))) {
    return(rep(NA_character_, nrow(evidence_rows)))
  }
  evidence_rows[[column_name]]
}

.select_assess_ecology_common_name <- function(evidence_rows, query_name, query_rank) {
  .select_assess_ecology_rank_metadata(
    evidence_rows = evidence_rows,
    query_name = query_name,
    query_rank = query_rank,
    column_name = "common_name",
    max_values = 1L
  )
}

.select_assess_ecology_rank_metadata <- function(
  evidence_rows,
  query_name,
  query_rank,
  column_name,
  max_values = 4L
) {
  if (!inherits(evidence_rows, "data.frame") ||
      nrow(evidence_rows) == 0 ||
      !(column_name %in% colnames(evidence_rows)) ||
      !("taxon_name" %in% colnames(evidence_rows)) ||
      !("taxon_rank" %in% colnames(evidence_rows))) {
    return(NA_character_)
  }

  query_key <- .taxon_match_key(query_name)
  query_rank <- tolower(trimws(as.character(query_rank)))
  if (!.is_non_empty_value(query_key) || !.is_non_empty_value(query_rank)) {
    return(NA_character_)
  }

  rows <- evidence_rows[
    .taxon_match_key(evidence_rows$taxon_name) == query_key &
      tolower(trimws(as.character(evidence_rows$taxon_rank))) == query_rank,
    ,
    drop = FALSE
  ]
  if (nrow(rows) == 0) {
    return(NA_character_)
  }

  .collapse_distinct_values(rows[[column_name]], max_values = max_values)
}

.assess_dominant_semicolon_value <- function(values, fallback = NULL) {
  tokens <- .assess_semicolon_terms(values)
  if (length(tokens) == 0 && !is.null(fallback)) {
    tokens <- .assess_semicolon_terms(fallback)
  }
  if (length(tokens) == 0) {
    return(list(value = NA_character_, prop = NA_real_))
  }
  counts <- as.data.frame(table(tokens), stringsAsFactors = FALSE)
  colnames(counts) <- c("value", "n")
  counts <- counts[order(-counts$n, counts$value), , drop = FALSE]
  list(
    value = as.character(counts$value[[1]]),
    prop = as.numeric(counts$n[[1]]) / sum(counts$n)
  )
}

.assess_semicolon_terms <- function(values) {
  values <- as.character(values)
  terms <- unlist(lapply(values, function(value) {
    split <- unlist(strsplit(value, "\\s*;\\s*"), use.names = FALSE)
    split <- trimws(split)
    unique(split[.is_non_empty_value(split)])
  }), use.names = FALSE)
  terms[.is_non_empty_value(terms)]
}

.assess_ecology_sum_numeric <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 0 || all(is.na(value))) {
    return(NA_real_)
  }
  sum(value, na.rm = TRUE)
}

.assess_ecology_min_numeric <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 0 || all(is.na(value))) {
    return(NA_real_)
  }
  min(value, na.rm = TRUE)
}

.assess_ecology_max_numeric <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 0 || all(is.na(value))) {
    return(NA_real_)
  }
  max(value, na.rm = TRUE)
}

.format_assess_ecology_summary <- function(summary_row) {
  taxa_returned <- if (is.null(summary_row) || nrow(summary_row) == 0) {
    NA_integer_
  } else {
    suppressWarnings(as.integer(summary_row$ecology_evidence_taxa_returned[[1]]))
  }
  if (is.na(taxa_returned) || taxa_returned == 0) {
    return(NA_character_)
  }
  paste0(
    "FishBase/SeaLifeBase ecology evidence: sources=",
    .format_summary_value(summary_row$ecology_evidence_sources[[1]]),
    "; taxa_returned=",
    taxa_returned,
    "; habitat_broad=",
    .format_summary_value(summary_row$ecology_evidence_habitat_broad[[1]]),
    "; depth_zone=",
    .format_summary_value(summary_row$ecology_evidence_depth_zone[[1]]),
    "; evidence_strength=",
    .format_summary_value(summary_row$ecology_evidence_strength[[1]]),
    "."
  )
}

.taxon_evidence_rank_summary <- function(rank, rank_values, taxon_names) {
  idx <- which(rank_values == rank & .is_non_empty_value(taxon_names))
  if (length(idx) == 0) {
    return(list(n = 0L, taxa = NA_character_))
  }
  names <- unique(taxon_names[idx][.is_non_empty_value(taxon_names[idx])])
  list(
    n = length(names),
    taxa = if (length(names) == 0) NA_character_ else paste(names, collapse = " | ")
  )
}

.infer_assess_taxon_evidence_region <- function(evidence_rows, region_label) {
  if (!is.null(region_label) &&
      length(region_label) > 0 &&
      any(.is_non_empty_value(region_label))) {
    return(as.character(region_label))
  }
  region <- .optional_evidence_column(evidence_rows, "region")
  .collapse_distinct_values(region, max_values = 8)
}

.add_assess_taxon_evidence_columns <- function(result, taxon_summary) {
  taxon_columns <- .assess_taxa_user_evidence_output_columns()
  for (column in taxon_columns) {
    if (!(column %in% colnames(result))) {
      result[[column]] <- NA
    }
  }
  if (is.null(taxon_summary) || nrow(taxon_summary) == 0 || !("query_id" %in% colnames(result))) {
    return(result)
  }

  lookup <- match(as.character(result$query_id), as.character(taxon_summary$query_id))
  matched <- !is.na(lookup)
  for (column in intersect(taxon_columns, colnames(taxon_summary))) {
    result[[column]][matched] <- taxon_summary[[column]][lookup[matched]]
  }
  result
}

summarise_obis_checklist <- function(chk, query_name, region = NA_character_) {
  query_name <- as.character(query_name)
  region <- as.character(region)
  if (length(region) == 0 || is.na(region) || !nzchar(trimws(region))) {
    region <- NA_character_
  }

  if (is.null(chk) || nrow(as.data.frame(chk)) == 0) {
    return(.empty_obis_checklist_summary(
      query_name = query_name,
      region = region,
      obis_region_status = "no_obis_records"
    ))
  }

  chk <- janitor::clean_names(as.data.frame(chk, stringsAsFactors = FALSE))
  records <- .obis_numeric_column(chk, "records")
  rank_values <- tolower(.obis_character_column(chk, c("taxon_rank", "rank")))
  scientific_names <- .obis_character_column(
    chk,
    c("scientific_name", "scientificname", "taxon_name")
  )

  out <- .empty_obis_checklist_summary(
    query_name = query_name,
    region = region,
    obis_region_status = "obis_records_found"
  )
  out$obis_taxa_returned <- nrow(chk)
  out$obis_total_records_including_descendants <- sum(records, na.rm = TRUE)
  out$obis_environment <- .obis_environment_summary(chk)

  for (rank in .assess_taxa_obis_ranks()) {
    rank_summary <- .obis_rank_summary(
      rank = rank,
      rank_values = rank_values,
      scientific_names = scientific_names,
      records = records
    )
    out[[paste0("obis_n_", rank)]] <- rank_summary$n
    out[[paste0("obis_taxa_", rank)]] <- rank_summary$taxa
  }

  out
}

fetch_obis_checklist_evidence <- function(
  query_taxa,
  obis_geometry,
  obis_region_label,
  obis_region_relation = c("direct", "supporting"),
  verbose = FALSE
) {
  obis_region_relation <- match.arg(obis_region_relation)
  .require_robis()
  if (!inherits(query_taxa, "data.frame") ||
      !all(c("query_name", "query_rank") %in% colnames(query_taxa))) {
    stop(
      "`query_taxa` must be a data.frame containing `query_name` and `query_rank`.",
      call. = FALSE
    )
  }
  if (nrow(query_taxa) == 0) {
    return(.empty_obis_checklist_evidence())
  }

  rows <- lapply(seq_len(nrow(query_taxa)), function(row_index) {
    query_name <- as.character(query_taxa$query_name[[row_index]])
    query_id <- if ("query_id" %in% colnames(query_taxa)) {
      as.character(query_taxa$query_id[[row_index]])
    } else {
      NA_character_
    }
    if (isTRUE(verbose)) {
      message("assess_taxa_evidence(): querying OBIS checklist for ", query_name, ".")
    }

    summary <- tryCatch(
      summarise_obis_checklist(
        chk = .obis_checklist(scientificname = query_name, geometry = obis_geometry),
        query_name = query_name,
        region = obis_region_label
      ),
      error = function(error) {
        out <- .empty_obis_checklist_summary(
          query_name = query_name,
          region = obis_region_label,
          obis_region_status = "obis_lookup_failed"
        )
        out$obis_error <- conditionMessage(error)
        out
      }
    )
    summary$query_id <- query_id
    summary$query_rank <- as.character(query_taxa$query_rank[[row_index]])
    summary$obis_region_relation <- obis_region_relation
    if (!("obis_error" %in% colnames(summary))) {
      summary$obis_error <- NA_character_
    }
    summary
  })

  out <- .bind_rows_fill(rows)
  rownames(out) <- NULL
  out
}

.require_robis <- function() {
  if (!requireNamespace("robis", quietly = TRUE)) {
    stop(
      "OBIS checklist evidence requires the suggested package `robis`. ",
      "Install it with install.packages(\"robis\") or omit `\"obis\"` from ",
      "`evidence_sources`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.obis_checklist <- function(scientificname, geometry) {
  robis::checklist(scientificname = scientificname, geometry = geometry)
}

.empty_obis_checklist_evidence <- function() {
  out <- .empty_obis_checklist_summary(character(), character())
  out$query_id <- character()
  out$query_rank <- character()
  out$obis_region_relation <- character()
  out$obis_error <- character()
  out
}

.empty_obis_checklist_summary <- function(
  query_name,
  region,
  obis_region_status = "no_obis_records"
) {
  out <- data.frame(
    query_name = as.character(query_name),
    region = as.character(region),
    obis_region_status = as.character(obis_region_status),
    obis_taxa_returned = 0L,
    obis_total_records_including_descendants = 0,
    obis_environment = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (rank in .assess_taxa_obis_ranks()) {
    out[[paste0("obis_n_", rank)]] <- 0L
    out[[paste0("obis_taxa_", rank)]] <- NA_character_
  }
  out
}

.assess_taxa_obis_ranks <- function() {
  c("phylum", "class", "order", "family", "genus", "species")
}

.assess_taxa_obis_output_columns <- function() {
  unlist(
    c(
      "obis_region_status",
      "obis_taxa_returned",
      "obis_total_records_including_descendants",
      "obis_environment",
      lapply(.assess_taxa_obis_ranks(), function(rank) {
        c(paste0("obis_n_", rank), paste0("obis_taxa_", rank))
      })
    ),
    use.names = FALSE
  )
}

.obis_numeric_column <- function(x, column_name) {
  if (!(column_name %in% colnames(x))) {
    return(rep(0, nrow(x)))
  }
  suppressWarnings(as.numeric(x[[column_name]]))
}

.obis_character_column <- function(x, column_names) {
  column <- intersect(column_names, colnames(x))
  if (length(column) == 0) {
    return(rep(NA_character_, nrow(x)))
  }
  as.character(x[[column[[1]]]])
}

.obis_environment_summary <- function(chk) {
  flags <- c(
    marine = "is_marine",
    brackish = "is_brackish",
    freshwater = "is_freshwater",
    terrestrial = "is_terrestrial"
  )
  values <- names(flags)[vapply(flags, function(column) {
    column %in% colnames(chk) && any(.obis_truthy(chk[[column]]), na.rm = TRUE)
  }, logical(1))]
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(values, collapse = "; ")
}

.obis_truthy <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value %in% c("true", "t", "1", "yes", "y")
}

.obis_rank_summary <- function(rank, rank_values, scientific_names, records) {
  idx <- which(rank_values == rank & .is_non_empty_value(scientific_names))
  if (length(idx) == 0) {
    return(list(n = 0L, taxa = NA_character_))
  }
  rank_records <- records[idx]
  rank_names <- scientific_names[idx]
  rank_names <- rank_names[order(rank_records, decreasing = TRUE, na.last = TRUE)]
  rank_names <- unique(rank_names[.is_non_empty_value(rank_names)])
  list(
    n = length(rank_names),
    taxa = if (length(rank_names) == 0) NA_character_ else paste(rank_names, collapse = " | ")
  )
}

.apply_assess_taxa_regional_evidence <- function(
  result,
  unique_query_table,
  user_taxon_evidence,
  obis_evidence,
  expected_environment,
  expected_habitat,
  expected_region,
  taxon_evidence_region_relation,
  taxon_evidence_region_label,
  obis_region_relation
) {
  if (!inherits(result, "data.frame") || nrow(result) == 0) {
    return(result)
  }
  taxon_summary <- summarise_assess_taxon_evidence_matches(
    unique_query_table = unique_query_table,
    taxon_evidence = user_taxon_evidence,
    taxon_evidence_region_label = taxon_evidence_region_label,
    taxon_evidence_region_relation = taxon_evidence_region_relation
  )
  if (!is.null(taxon_summary) && nrow(taxon_summary) > 0) {
    result <- .add_assess_taxon_evidence_columns(result, taxon_summary)
  }
  if (is.null(expected_region)) {
    result <- .apply_assess_obis_environment_fallback(
      result = result,
      obis_evidence = obis_evidence,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat
    )
    result <- .update_assess_action_after_region_status(result)
    if (!is.null(obis_evidence) && nrow(obis_evidence) > 0) {
      result <- .add_assess_obis_columns(result, obis_evidence)
    }
    return(result)
  }

  taxon_region <- .assess_taxon_evidence_region_support(taxon_summary)
  obis_region <- .assess_obis_region_support(
    unique_query_table = unique_query_table,
    obis_evidence = obis_evidence,
    obis_region_relation = obis_region_relation
  )
  support <- merge(
    taxon_region,
    obis_region,
    by = "query_id",
    all = TRUE,
    sort = FALSE
  )

  support$taxon_evidence_region_status[is.na(support$taxon_evidence_region_status)] <- "no_distribution_evidence"
  support$obis_expected_region_status[is.na(support$obis_expected_region_status)] <- "no_distribution_evidence"
  support$combined_region_status <- mapply(
    .highest_assess_region_status,
    support$taxon_evidence_region_status,
    support$obis_expected_region_status,
    MoreArgs = list(current_status = "no_distribution_evidence"),
    USE.NAMES = FALSE
  )

  support_lookup <- stats::setNames(support$combined_region_status, support$query_id)
  original_status <- as.character(result$expected_region_status)
  result$expected_region_status <- vapply(seq_len(nrow(result)), function(row_index) {
    query_id <- as.character(result$query_id[[row_index]])
    status <- .normalise_assess_region_status(original_status[[row_index]])
    if (query_id %in% names(support_lookup)) {
      status <- .highest_assess_region_status(status, support_lookup[[query_id]])
    }
    status
  }, character(1))

  result <- .apply_assess_obis_environment_fallback(
    result = result,
    obis_evidence = obis_evidence,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat
  )
  result <- .update_assess_action_after_region_status(result)
  if (!is.null(obis_evidence) && nrow(obis_evidence) > 0) {
    result <- .add_assess_obis_columns(result, obis_evidence)
    result <- .append_assess_obis_provenance(result, obis_evidence)
  }
  result
}

.apply_assess_ecology_evidence <- function(
  result,
  unique_query_table,
  ecology_evidence,
  expected_environment,
  expected_habitat,
  expected_habitat_type,
  expected_region
) {
  if (!inherits(result, "data.frame") || nrow(result) == 0 || is.null(ecology_evidence)) {
    return(result)
  }

  ecology_summary <- summarise_assess_ecology_evidence_matches(
    unique_query_table = unique_query_table,
    ecology_evidence = ecology_evidence
  )
  if (is.null(ecology_summary) || nrow(ecology_summary) == 0) {
    return(result)
  }

  original_environment_status <- if ("expected_environment_status" %in% colnames(result)) {
    as.character(result$expected_environment_status)
  } else {
    NULL
  }
  original_region_status <- if ("expected_region_status" %in% colnames(result)) {
    as.character(result$expected_region_status)
  } else {
    NULL
  }

  result <- .add_assess_ecology_columns(result, ecology_summary)
  if (!("query_id" %in% colnames(result))) {
    return(result)
  }

  lookup <- match(as.character(result$query_id), as.character(ecology_summary$query_id))
  matched <- which(!is.na(lookup) & suppressWarnings(as.integer(ecology_summary$ecology_evidence_taxa_returned[lookup])) > 0)
  if (length(matched) == 0) {
    return(result)
  }

  for (row_index in matched) {
    ecology_row <- ecology_summary[lookup[[row_index]], , drop = FALSE]
    result <- .apply_assess_ecology_row(
      result = result,
      row_index = row_index,
      ecology_row = ecology_row,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_habitat_type = expected_habitat_type,
      expected_region = expected_region
    )
  }

  if (!is.null(original_environment_status)) {
    result$expected_environment_status <- original_environment_status
  }
  if (!is.null(original_region_status)) {
    result$expected_region_status <- original_region_status
  }
  result <- .update_assess_action_after_region_status(result)
  if (!is.null(original_environment_status)) {
    result$expected_environment_status <- original_environment_status
  }
  if (!is.null(original_region_status)) {
    result$expected_region_status <- original_region_status
  }
  .refresh_assess_rationale_final_status_text(result)
}

.add_assess_ecology_columns <- function(result, ecology_summary) {
  ecology_columns <- .assess_taxa_ecology_output_columns()
  for (column in ecology_columns) {
    if (!(column %in% colnames(result))) {
      result[[column]] <- NA
    }
  }
  if (is.null(ecology_summary) || nrow(ecology_summary) == 0 || !("query_id" %in% colnames(result))) {
    return(result)
  }

  lookup <- match(as.character(result$query_id), as.character(ecology_summary$query_id))
  matched <- !is.na(lookup)
  for (column in intersect(ecology_columns, colnames(ecology_summary))) {
    result[[column]][matched] <- ecology_summary[[column]][lookup[matched]]
  }
  result
}

.apply_assess_ecology_row <- function(
  result,
  row_index,
  ecology_row,
  expected_environment,
  expected_habitat,
  expected_habitat_type,
  expected_region
) {
  current_action <- as.character(result$recommended_action[[row_index]])
  habitat_status <- .assess_ecology_habitat_status(
    ecology_row = ecology_row,
    expected_habitat = expected_habitat,
    expected_habitat_type = expected_habitat_type,
    current_status = as.character(result$expected_habitat_status[[row_index]])
  )
  if (!identical(habitat_status, as.character(result$expected_habitat_status[[row_index]]))) {
    result$expected_habitat_status[[row_index]] <- habitat_status
  }

  result$evidence_summary[[row_index]] <- .append_assess_text(
    result$evidence_summary[[row_index]],
    as.character(ecology_row$ecology_evidence_summary[[1]])
  )
  result$rationale[[row_index]] <- .append_assess_text(
    result$rationale[[row_index]],
    .format_assess_ecology_rationale(ecology_row, habitat_status)
  )
  result$references[[row_index]] <- .append_assess_text(
    result$references[[row_index]],
    "FishBase/SeaLifeBase local ecology DB",
    separator = "; "
  )
  result$evidence_sources[[row_index]] <- .append_assess_source(
    result$evidence_sources[[row_index]],
    as.character(ecology_row$ecology_evidence_sources[[1]])
  )
  if ("evidence_basis" %in% colnames(result) &&
      identical(as.character(result$evidence_basis[[row_index]]), "conservative_reasoning_only")) {
    result$evidence_basis[[row_index]] <- "local_evidence"
  }

  if (current_action %in% c("exclude", "flag_possible_exclusion") &&
      !identical(habitat_status, "possible_unlikely")) {
    return(result)
  }

  interpretation <- .assess_ecology_interpretation(
    env_status = as.character(result$expected_environment_status[[row_index]]),
    habitat_status = habitat_status,
    region_status = as.character(result$expected_region_status[[row_index]]),
    current_action = current_action
  )
  if (.assess_ecology_pelagic_deep_benthic_signal(ecology_row, expected_habitat_type, habitat_status)) {
    interpretation$ecological_status <- "plausible"
    interpretation$occurrence_interpretation <- "possible_transient_or_allochthonous"
    if (interpretation$recommended_action %in% c("flag_for_review", "retain")) {
      interpretation$recommended_action <- "retain"
    }
    result$rationale[[row_index]] <- .append_assess_text(
      result$rationale[[row_index]],
      "FishBase/SeaLifeBase ecology indicates pelagic ecology in a deep-sea benthic context. The taxon is retained as plausible biological signal when other evidence supports marine/regional plausibility, but interpreted as possible transient or allochthonous material rather than a resident benthic organism."
    )
  }
  result$ecological_status[[row_index]] <- interpretation$ecological_status
  result$occurrence_interpretation[[row_index]] <- interpretation$occurrence_interpretation
  result$recommended_action[[row_index]] <- interpretation$recommended_action
  result
}

.assess_ecology_pelagic_deep_benthic_signal <- function(ecology_row, expected_habitat_type, habitat_status) {
  if (!identical(expected_habitat_type, "deep_sea_benthic") ||
      !identical(as.character(habitat_status), "possible_likely")) {
    return(FALSE)
  }
  habitat_values <- .assess_ecology_split_values(
    c(
      ecology_row$ecology_evidence_habitat_broad,
      ecology_row$ecology_evidence_dominant_habitat_broad
    )
  )
  length(habitat_values) > 0 && all(habitat_values %in% "pelagic")
}

.assess_ecology_habitat_status <- function(
  ecology_row,
  expected_habitat,
  expected_habitat_type,
  current_status
) {
  current_status <- as.character(current_status)
  if (is.null(expected_habitat) || !.is_non_empty_value(expected_habitat)) {
    return(current_status)
  }
  if (current_status %in% c("compatible", "incompatible")) {
    return(current_status)
  }

  context <- .assess_ecology_expected_habitat_context(expected_habitat)
  habitat_values <- .assess_ecology_split_values(
    c(
      ecology_row$ecology_evidence_habitat_broad,
      ecology_row$ecology_evidence_dominant_habitat_broad,
      ecology_row$ecology_evidence_position
    )
  )
  depth_values <- .assess_ecology_split_values(ecology_row$ecology_evidence_depth_zone)
  text_values <- tolower(paste(c(
    ecology_row$ecology_evidence_rank_habitat_note,
    ecology_row$ecology_evidence_summary
  ), collapse = " "))

  if (.flag_taxa_values_contain_context(c(habitat_values, depth_values, text_values), expected_habitat)) {
    return("compatible")
  }
  if ("deep" %in% context && .assess_ecology_depth_supports_deep(depth_values)) {
    return("possible_likely")
  }
  if ("benthic" %in% context && any(habitat_values %in% c("benthic", "demersal", "benthopelagic"))) {
    return("compatible")
  }
  if ("pelagic" %in% context && any(habitat_values %in% c("pelagic", "benthopelagic"))) {
    return("compatible")
  }
  if ("host" %in% context && any(habitat_values %in% c("host_associated"))) {
    return("compatible")
  }
  if (identical(expected_habitat_type, "deep_sea_benthic") &&
      length(habitat_values) > 0 &&
      all(habitat_values %in% "pelagic")) {
    return("possible_likely")
  }
  if ("deep" %in% context && any(depth_values %in% c("shallow_only"))) {
    return("possible_unlikely")
  }
  if (any(c("deep", "benthic") %in% context) &&
      length(habitat_values) > 0 &&
      all(habitat_values %in% c("pelagic", "host_associated"))) {
    return("possible_unlikely")
  }
  if ("pelagic" %in% context &&
      length(habitat_values) > 0 &&
      all(habitat_values %in% c("benthic", "demersal", "host_associated"))) {
    return("possible_unlikely")
  }

  current_status
}

.assess_ecology_expected_habitat_context <- function(expected_habitat) {
  habitat <- tolower(as.character(expected_habitat))
  context <- character()
  if (grepl("deep sea|deep-sea|abyss|hadal|bathyal|aphotic", habitat)) {
    context <- c(context, "deep")
  }
  if (grepl("benthic|benthos|sediment|seafloor|bottom|demersal", habitat)) {
    context <- c(context, "benthic")
  }
  if (grepl("pelagic|water column|plankton|nekton", habitat)) {
    context <- c(context, "pelagic")
  }
  if (grepl("host|parasite|symbiont|epiphyt", habitat)) {
    context <- c(context, "host")
  }
  unique(context)
}

.assess_ecology_depth_supports_deep <- function(depth_values) {
  any(depth_values %in% c(
    "mesophotic_or_upper_bathyal_possible",
    "bathyal_possible",
    "abyssal_possible",
    "hadal_or_abyssal_possible"
  ))
}

.assess_ecology_split_values <- function(value) {
  value <- as.character(value)
  value <- value[.is_non_empty_value(value)]
  if (length(value) == 0) {
    return(character())
  }
  values <- unlist(strsplit(value, "\\s*;\\s*|\\s*,\\s*|\\s*/\\s*"), use.names = FALSE)
  values <- tolower(trimws(values))
  unique(values[.is_non_empty_value(values) & values != "unknown"])
}

.format_assess_ecology_rationale <- function(ecology_row, habitat_status) {
  paste0(
    "FishBase/SeaLifeBase ecology evidence was matched for this taxon and is used only as habitat/ecology metadata. ",
    "It does not modify expected_environment_status or expected_region_status. ",
    "Ecology habitat interpretation from this layer: ",
    habitat_status,
    "."
  )
}

.assess_ecology_interpretation <- function(
  env_status,
  habitat_status,
  region_status,
  current_action
) {
  if (identical(habitat_status, "possible_unlikely")) {
    return(list(
      ecological_status = "unlikely_resident",
      occurrence_interpretation = "possible_transient_or_allochthonous",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  if (
    env_status %in% c("compatible", "mixed_within_rank", "possible_likely") &&
      habitat_status %in% c("compatible", "possible_likely", "unknown") &&
      region_status %in% c("known_in_region", "known_near_region", "plausible_in_region")
  ) {
    fully_known <- identical(habitat_status, "compatible") &&
      region_status %in% c("known_in_region", "known_near_region")
    return(list(
      ecological_status = if (fully_known) "compatible" else "plausible",
      occurrence_interpretation = if (fully_known) "expected_resident" else "plausible_resident",
      recommended_action = "retain"
    ))
  }

  if (
    env_status %in% c("compatible", "mixed_within_rank", "possible_likely") &&
      habitat_status %in% c("compatible", "possible_likely") &&
      region_status %in% c("not_assessed", "unknown")
  ) {
    return(list(
      ecological_status = if (identical(habitat_status, "compatible")) "compatible" else "plausible",
      occurrence_interpretation = if (identical(habitat_status, "compatible")) "expected_resident" else "plausible_resident",
      recommended_action = "retain"
    ))
  }

  list(
    ecological_status = if (habitat_status %in% c("compatible", "possible_likely")) "plausible" else "unknown",
    occurrence_interpretation = if (habitat_status %in% c("compatible", "possible_likely")) "plausible_resident" else "uncertain",
    recommended_action = if (current_action %in% c("retain", "flag_for_review")) current_action else "flag_for_review"
  )
}

.apply_assess_deep_sea_benthic_vertebrate_rule <- function(result, expected_habitat_type) {
  if (!identical(expected_habitat_type, "deep_sea_benthic") ||
      !inherits(result, "data.frame") ||
      nrow(result) == 0) {
    return(result)
  }
  required <- c(
    "taxon_name",
    "taxon_rank",
    "lineage",
    "expected_environment_status",
    "expected_habitat_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale"
  )
  if (!all(required %in% colnames(result))) {
    return(result)
  }

  target <- .assess_deep_sea_benthic_air_breathing_class(result) %in% c("aves", "reptilia") &
    !(result$recommended_action %in% "exclude") &
    !(result$expected_environment_status %in% "incompatible")
  if (!any(target, na.rm = TRUE)) {
    return(result)
  }

  result$expected_habitat_status[target] <- "possible_unlikely"
  result$ecological_status[target] <- "habitat_mismatch_or_allochthonous"
  result$occurrence_interpretation[target] <- "possible_transient_or_allochthonous"
  result$recommended_action[target] <- "flag_possible_exclusion"
  result$rationale[target] <- vapply(
    result$rationale[target],
    .append_assess_text,
    character(1),
    addition = "The expected habitat type is deep_sea_benthic and the taxon lineage indicates Aves or Reptilia. This is treated as a resident-habitat mismatch and possible transient/allochthonous signal, not as a blanket environment exclusion. Mammalia are not included in this blanket rule because several marine mammals dive deeply."
  )
  .refresh_assess_rationale_final_status_text(result)
}

.assess_deep_sea_benthic_air_breathing_class <- function(result) {
  text <- paste(
    as.character(result$taxon_name),
    as.character(result$taxon_rank),
    as.character(result$lineage),
    sep = " "
  )
  text <- tolower(text)
  out <- rep(NA_character_, length(text))
  out[grepl("class=aves|\\baves\\b|\\bbird\\b", text)] <- "aves"
  out[grepl("class=reptilia|\\breptilia\\b|\\breptile\\b", text)] <- "reptilia"
  out[grepl("class=mammalia|\\bmammalia\\b|\\bmammal\\b", text)] <- "mammalia"
  out
}

.assess_taxon_evidence_region_support <- function(taxon_summary) {
  if (is.null(taxon_summary) || nrow(taxon_summary) == 0) {
    return(data.frame(
      query_id = character(),
      taxon_evidence_region_status = character(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  relation <- as.character(taxon_summary$taxon_evidence_region_relation)
  has_evidence <- suppressWarnings(as.integer(taxon_summary$taxon_evidence_taxa_returned)) > 0
  has_region_context <- .is_non_empty_value(taxon_summary$taxon_evidence_region)
  status <- rep("no_distribution_evidence", nrow(taxon_summary))
  status[has_evidence & has_region_context & relation == "direct"] <- "known_in_region"
  status[has_evidence & has_region_context & relation == "supporting"] <- "plausible_in_region"
  data.frame(
    query_id = as.character(taxon_summary$query_id),
    taxon_evidence_region_status = status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.assess_obis_region_support <- function(
  unique_query_table,
  obis_evidence,
  obis_region_relation
) {
  out <- data.frame(
    query_id = as.character(unique_query_table$query_id),
    obis_expected_region_status = rep("no_distribution_evidence", nrow(unique_query_table)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (is.null(obis_evidence) || nrow(obis_evidence) == 0) {
    return(out)
  }
  support_status <- if (identical(obis_region_relation, "direct")) {
    "known_in_region"
  } else {
    "plausible_in_region"
  }
  rows_with_records <- as.character(obis_evidence$query_id)[
    .obis_summary_has_records(obis_evidence)
  ]
  out$obis_expected_region_status[out$query_id %in% rows_with_records] <- support_status
  out
}

.apply_assess_obis_environment_fallback <- function(
  result,
  obis_evidence,
  expected_environment,
  expected_habitat
) {
  if (is.null(obis_evidence) || nrow(obis_evidence) == 0 ||
      !all(c("query_id", "expected_environment_status") %in% colnames(result)) ||
      !("obis_environment" %in% colnames(obis_evidence))) {
    return(result)
  }

  lookup <- match(as.character(result$query_id), as.character(obis_evidence$query_id))
  rows <- which(!is.na(lookup) &
    result$expected_environment_status == "unknown" &
    .obis_summary_has_records(obis_evidence[lookup, , drop = FALSE]))
  if (length(rows) == 0) {
    return(result)
  }

  for (row_index in rows) {
    obis_row <- obis_evidence[lookup[[row_index]], , drop = FALSE]
    env_values <- .flag_taxa_evidence_values(
      data.frame(environment = as.character(obis_row$obis_environment[[1]]), stringsAsFactors = FALSE),
      "environment"
    )
    supports_environment <- .flag_taxa_values_contain_context(env_values, expected_environment)
    if (!isTRUE(supports_environment)) {
      next
    }

    query_rank <- if ("taxon_rank" %in% colnames(result)) {
      as.character(result$taxon_rank[[row_index]])
    } else if ("query_rank" %in% colnames(result)) {
      as.character(result$query_rank[[row_index]])
    } else {
      NA_character_
    }
    is_species <- identical(tolower(query_rank), "species")
    env_status <- .flag_taxa_prelim_environment_status(
      env_values = env_values,
      expected_environment = expected_environment,
      supports_environment = supports_environment,
      contradicts_environment = FALSE,
      is_species = is_species
    )
    result$expected_environment_status[[row_index]] <- env_status

    evidence_rows <- data.frame(
      environment = paste(env_values, collapse = "; "),
      evidence_summary = "OBIS checklist environment flags used as fallback because no matched user/WoRMS environment evidence was available.",
      reference = "OBIS",
      stringsAsFactors = FALSE
    )
    interpretation <- .flag_taxa_prelim_interpretation(
      env_status = env_status,
      habitat_status = as.character(result$expected_habitat_status[[row_index]]),
      region_status = as.character(result$expected_region_status[[row_index]]),
      query_rank = query_rank,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      evidence_rows = evidence_rows,
      is_species = is_species,
      has_positive_environment_incompatibility = FALSE
    )
    result$ecological_status[[row_index]] <- interpretation$ecological_status
    result$occurrence_interpretation[[row_index]] <- interpretation$occurrence_interpretation
    result$recommended_action[[row_index]] <- interpretation$recommended_action
    result$rationale[[row_index]] <- .append_assess_text(
      result$rationale[[row_index]],
      paste0(
        "OBIS environment flags (",
        paste(env_values, collapse = "; "),
        ") were used as fallback environment evidence because no matched user/WoRMS environment evidence was available."
      )
    )
    result$evidence_summary[[row_index]] <- .append_assess_text(
      result$evidence_summary[[row_index]],
      paste0(
        "OBIS environment fallback: environment=",
        paste(env_values, collapse = "; "),
        "."
      )
    )
    result$evidence_sources[[row_index]] <- .append_assess_source(
      result$evidence_sources[[row_index]],
      "obis"
    )
    if ("evidence_basis" %in% colnames(result) &&
        identical(result$evidence_basis[[row_index]], "conservative_reasoning_only")) {
      result$evidence_basis[[row_index]] <- "local_evidence"
    }
  }

  result
}

.normalise_assess_region_status <- function(status) {
  status <- as.character(status)
  if (identical(status, "known_in_region")) {
    return("known_in_region")
  }
  if (status %in% c("plausible_in_region", "known_near_region")) {
    return("plausible_in_region")
  }
  "no_distribution_evidence"
}

.highest_assess_region_status <- function(..., current_status = NULL) {
  statuses <- c(current_status, unlist(list(...), use.names = FALSE))
  statuses <- vapply(statuses, .normalise_assess_region_status, character(1))
  priority <- c(
    no_distribution_evidence = 1L,
    plausible_in_region = 2L,
    known_in_region = 3L
  )
  statuses[[which.max(priority[statuses])]]
}

.update_assess_action_after_region_status <- function(result) {
  if (!all(c(
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "recommended_action",
    "ecological_status",
    "occurrence_interpretation"
  ) %in% colnames(result))) {
    return(result)
  }

  retain <- result$recommended_action %in% c("flag_for_review", "retain") &
    result$expected_habitat_status %in% c("compatible", "possible_likely", "unknown") &
    result$expected_region_status %in% c("known_in_region", "known_near_region", "plausible_in_region") &
    result$expected_environment_status %in% c("compatible", "mixed_within_rank", "possible_likely")
  if (any(retain, na.rm = TRUE)) {
    preserve_transient <- retain &
      result$occurrence_interpretation %in% "possible_transient_or_allochthonous"
    result$recommended_action[retain] <- "retain"
    fully_known <- retain &
      result$expected_habitat_status == "compatible" &
      result$expected_region_status %in% c("known_in_region", "known_near_region")
    result$ecological_status[retain] <- "plausible"
    result$occurrence_interpretation[retain] <- "plausible_resident"
    result$ecological_status[fully_known] <- "compatible"
    result$occurrence_interpretation[fully_known] <- "expected_resident"
    result$ecological_status[preserve_transient] <- "plausible"
    result$occurrence_interpretation[preserve_transient] <- "possible_transient_or_allochthonous"

    unknown_habitat <- retain & result$expected_habitat_status == "unknown"
    result$rationale[unknown_habitat] <- vapply(
      result$rationale[unknown_habitat],
      .append_assess_text,
      character(1),
      addition = "No direct habitat-specific evidence was available in the deterministic evidence sources, so habitat is marked unknown. Because the taxon is environmentally compatible and known/plausible in the region, and no positive incompatibility signal was found, the taxon is retained."
    )

    broad_rank <- retain & tolower(as.character(result$taxon_rank)) %in% c("phylum", "class", "order", "family", "genus")
    result$rationale[broad_rank] <- vapply(
      result$rationale[broad_rank],
      .append_assess_text,
      character(1),
      addition = "The taxon is resolved at a broad rank, but available evidence supports environmental compatibility and regional plausibility. Broad rank alone is not treated as a reason for review."
    )

    supporting_region <- retain & result$expected_region_status == "plausible_in_region"
    result$rationale[supporting_region] <- vapply(
      result$rationale[supporting_region],
      .append_assess_text,
      character(1),
      addition = "Regional evidence is supporting rather than direct, so expected_region_status is plausible_in_region rather than known_in_region."
    )
  }

  .refresh_assess_rationale_final_status_text(result)
}

.refresh_assess_rationale_final_status_text <- function(result) {
  if (!inherits(result, "data.frame") || nrow(result) == 0 ||
      !("rationale" %in% colnames(result))) {
    return(result)
  }

  for (row_index in seq_len(nrow(result))) {
    rationale <- as.character(result$rationale[[row_index]])
    if (!.is_non_empty_value(rationale)) {
      next
    }
    if ("expected_environment_status" %in% colnames(result)) {
      rationale <- sub(
        "Environment status: [^.]+\\.",
        paste0("Environment status: ", result$expected_environment_status[[row_index]], "."),
        rationale
      )
    }
    if ("expected_habitat_status" %in% colnames(result)) {
      rationale <- sub(
        "Habitat status: [^.]+\\.",
        paste0("Habitat status: ", result$expected_habitat_status[[row_index]], "."),
        rationale
      )
    }
    if ("expected_region_status" %in% colnames(result)) {
      rationale <- sub(
        "Region status: [^.]+\\.",
        paste0("Region status: ", result$expected_region_status[[row_index]], "."),
        rationale
      )
    }
    if ("recommended_action" %in% colnames(result) &&
        identical(as.character(result$recommended_action[[row_index]]), "retain")) {
      rationale <- gsub(
        " The query rank is broad, but broad rank alone is not treated as a reason for review; this row is reviewed because evidence remains missing, conflicting, or insufficient for a plausible deterministic call.",
        "",
        rationale,
        fixed = TRUE
      )
    }
    result$rationale[[row_index]] <- rationale
  }

  result
}

.add_assess_obis_columns <- function(result, obis_evidence) {
  obis_columns <- .assess_taxa_obis_output_columns()
  for (column in obis_columns) {
    if (!(column %in% colnames(result))) {
      result[[column]] <- NA
    }
  }
  if (is.null(obis_evidence) || nrow(obis_evidence) == 0 || !("query_id" %in% colnames(result))) {
    return(result)
  }

  lookup <- match(as.character(result$query_id), as.character(obis_evidence$query_id))
  matched <- !is.na(lookup)
  for (column in intersect(obis_columns, colnames(obis_evidence))) {
    result[[column]][matched] <- obis_evidence[[column]][lookup[matched]]
  }
  result
}

.append_assess_obis_provenance <- function(result, obis_evidence) {
  if (is.null(obis_evidence) || nrow(obis_evidence) == 0 || !("query_id" %in% colnames(result))) {
    return(result)
  }

  lookup <- match(as.character(result$query_id), as.character(obis_evidence$query_id))
  matched <- !is.na(lookup)
  for (row_index in which(matched)) {
    obis_row <- obis_evidence[lookup[[row_index]], , drop = FALSE]
    summary <- .format_assess_obis_reference(obis_row)
    result$evidence_summary[[row_index]] <- .append_assess_text(
      result$evidence_summary[[row_index]],
      summary
    )
    result$rationale[[row_index]] <- .append_assess_text(
      result$rationale[[row_index]],
      .format_assess_obis_rationale(obis_row)
    )
    if (.obis_summary_has_records(obis_row)) {
      result$references[[row_index]] <- .append_assess_text(
        result$references[[row_index]],
        summary,
        separator = "; "
      )
      result$evidence_sources[[row_index]] <- .append_assess_source(
        result$evidence_sources[[row_index]],
        "obis"
      )
      if (identical(result$evidence_basis[[row_index]], "conservative_reasoning_only")) {
        result$evidence_basis[[row_index]] <- "local_evidence"
      }
    }
  }

  result
}

.obis_summary_has_records <- function(obis_row) {
  taxa_returned <- suppressWarnings(as.integer(obis_row$obis_taxa_returned))
  status <- as.character(obis_row$obis_region_status)
  !is.na(taxa_returned) & taxa_returned > 0 & status == "obis_records_found"
}

.format_assess_obis_reference <- function(obis_row) {
  query_name <- as.character(obis_row$query_name[[1]])
  region <- as.character(obis_row$region[[1]])
  records <- as.numeric(obis_row$obis_total_records_including_descendants[[1]])
  if (is.na(records)) {
    records <- 0
  }
  if (identical(as.character(obis_row$obis_region_status[[1]]), "obis_lookup_failed")) {
    error_message <- as.character(obis_row$obis_error[[1]])
    return(paste0(
      "OBIS checklist query for ",
      query_name,
      " in ",
      region,
      " failed: ",
      error_message
    ))
  }
  if (.obis_summary_has_records(obis_row)) {
    return(paste0(
      "OBIS checklist query for ",
      query_name,
      " in ",
      region,
      "; ",
      records,
      " records reported across returned taxa/descendants."
    ))
  }
  paste0(
    "OBIS checklist query for ",
    query_name,
    " in ",
    region,
    "; no OBIS checklist records were returned for the supplied geometry."
  )
}

.format_assess_obis_rationale <- function(obis_row) {
  if (.obis_summary_has_records(obis_row)) {
    return(paste0(
      .format_assess_obis_reference(obis_row),
      " OBIS is treated as regional occurrence support, not as proof of resident habitat use."
    ))
  }
  paste0(
    .format_assess_obis_reference(obis_row),
    " No OBIS checklist records were returned for the supplied geometry. This is treated as missing regional evidence, not evidence of absence."
  )
}

.append_assess_text <- function(existing, addition, separator = " ") {
  existing <- as.character(existing)
  addition <- as.character(addition)
  if (!.is_non_empty_value(addition)) {
    return(existing)
  }
  if (!.is_non_empty_value(existing)) {
    return(addition)
  }
  paste(existing, addition, sep = separator)
}

.append_assess_source <- function(existing, source) {
  existing <- as.character(existing)
  source <- as.character(source)
  if (!.is_non_empty_value(existing) || identical(existing, "none")) {
    return(source)
  }
  sources <- unique(c(unlist(strsplit(existing, "\\s*;\\s*")), source))
  paste(sources[.is_non_empty_value(sources)], collapse = "; ")
}
