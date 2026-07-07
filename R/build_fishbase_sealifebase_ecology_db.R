#' Download FishBase and SeaLifeBase parquet files
#'
#' @description
#' `download_fishbase_sealifebase_parquet()` downloads the FishBase and
#' SeaLifeBase species, ecology, and native taxonomy parquet files used by
#' `build_fishbase_sealifebase_ecology_db()`. The files are saved under
#' `output_dir` and are intended as local, uncommitted build inputs.
#'
#' @param output_dir Directory where the `fb/` and `slb/` parquet folders should
#'   be created. Defaults to `"data-raw/fishbase_parquet"`.
#' @param version FishBase/SeaLifeBase parquet release version. Defaults to
#'   `"25.04"`.
#' @param overwrite Logical scalar. If `FALSE`, existing files are not
#'   downloaded again.
#'
#' @return A data frame with source, table, URL, path, and downloaded status.
#' @export
#'
#' @examples
#' \dontrun{
#' download_fishbase_sealifebase_parquet()
#' }
download_fishbase_sealifebase_parquet <- function(
  output_dir = "data-raw/fishbase_parquet",
  version = "25.04",
  overwrite = FALSE
) {
  output_dir <- .fb_required_scalar_character(output_dir, "output_dir")
  version <- .fb_required_scalar_character(version, "version")
  overwrite <- .fb_logical_scalar(overwrite, "overwrite")
  .require_curl_for_fishbase_download()

  manifest <- .fishbase_sealifebase_parquet_manifest(output_dir, version)
  for (row_index in seq_len(nrow(manifest))) {
    path <- manifest$path[[row_index]]
    if (file.exists(path) && !isTRUE(overwrite)) {
      manifest$downloaded[[row_index]] <- FALSE
      next
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      curl::curl_download(
        url = manifest$url[[row_index]],
        destfile = path,
        quiet = FALSE
      ),
      error = function(error) {
        stop(
          "Failed to download FishBase/SeaLifeBase parquet file from:\n",
          manifest$url[[row_index]],
          "\nDestination:\n",
          path,
          "\nOriginal error: ",
          conditionMessage(error),
          "\nDownload the file manually, for example with command-line curl:\n",
          "curl -L -o ",
          shQuote(path),
          " ",
          shQuote(manifest$url[[row_index]]),
          call. = FALSE
        )
      }
    )
    manifest$downloaded[[row_index]] <- TRUE
  }

  manifest
}

#' Build a FishBase/SeaLifeBase ecology evidence database
#'
#' @description
#' `build_fishbase_sealifebase_ecology_db()` reads local FishBase and
#' SeaLifeBase parquet exports, reconstructs native FishBase/SeaLifeBase
#' taxonomy locally, cleans species and ecology tables separately by source,
#' collapses ecology flags into compact metadata columns, joins by `source` and
#' `spec_code`, and returns species-level and higher-rank ecology evidence
#' tables. Species-level rows are uniquely identified by `source + spec_code`.
#'
#' This function does not use `rfishbase` and does not make live API calls. If
#' parquet files are missing and `download_if_missing = TRUE`, it downloads the
#' current Source Cooperative parquet files with `curl::curl_download()`.
#'
#' FishBase/SeaLifeBase native taxonomy is used rather than NCBI taxonomy.
#' FishBase kingdom/phylum are inferred as `Metazoa`/`Chordata` because the
#' native FishBase taxonomy tables are fish-focused and do not carry those ranks
#' explicitly. Only standard ranks (`kingdom`, `phylum`, `class`, `order`,
#' `family`, `genus`, and `species`) are emitted as `taxon_rank` rows.
#'
#' Environment evidence is represented in a WoRMS-like collapsed form such as
#' `"marine; brackish; freshwater"`. `ecology_habitat_broad` is derived only
#' from the species-table `DemersPelag` field. Ecology-table flags are retained
#' as compact metadata columns such as `ecology_zone`,
#' `ecology_water_column_zone`, `ecology_mobility`, `ecology_size_class`,
#' `ecology_substrate`, and `ecology_special_habitat`, but they do not determine
#' `ecology_habitat_broad`.
#'
#' `depth_min` and `depth_max` come from the full reported depth range;
#' `common_depth_min` and `common_depth_max` come from the common depth range
#' when available. `ecology_depth_zone` is based on `depth_max`, not
#' `common_depth_max`, because it is intended to capture possible deep
#' occurrence rather than common occurrence.
#'
#' Aggregated semicolon-separated ecology terms are split, trimmed, and
#' de-duplicated before being returned. Dominant ecology fields such as
#' `dominant_environment`, `dominant_habitat_broad`, and
#' `dominant_depth_zone` contain one deterministic value only, with the paired
#' `*_prop` column giving the dominant term count divided by all usable split
#' terms.
#'
#' Higher-rank rows are aggregated from species-level ecology evidence, while
#' native rank-table metadata is retained in fields such as `common_name`,
#' `rank_water_salinity`, `rank_habitat_note`, `rank_distribution`, and
#' `rank_etymology`.
#' `rank_habitat_note` is native rank-table text metadata and should not be
#' confused with aggregated `ecology_habitat_broad`. For higher-rank rows,
#' `ecology_evidence_strength` summarises support for the dominant broad habitat:
#' `high` means at least 10 species with habitat evidence and a dominant
#' proportion of at least 0.8; `moderate` means at least 5 species and a
#' dominant proportion of at least 0.6; `mixed` means multiple broad habitats
#' with no moderate/high dominance; `low` means fewer than 5 species with some
#' habitat evidence; and `unknown` means no useful species-level broad habitat
#' evidence.
#'
#' The returned list includes diagnostics with source-level row counts and
#' unmatched ecology/taxonomy `spec_code` counts. The final species-level row
#' count is checked against the combined FishBase and SeaLifeBase species row
#' count. FishBase kingdom/phylum are inferred as `Metazoa`/`Chordata` because
#' the native FishBase taxonomy tables are fish-focused and do not carry those
#' ranks explicitly.
#'
#' @param fishbase_species_path Path to the FishBase `species.parquet` file.
#' @param fishbase_ecology_path Path to the FishBase `ecology.parquet` file.
#' @param fishbase_genera_path Path to the FishBase `genera.parquet` file.
#' @param fishbase_families_path Path to the FishBase `families.parquet` file.
#' @param fishbase_orders_path Path to the FishBase `orders.parquet` file.
#' @param fishbase_classes_path Path to the FishBase `classes.parquet` file.
#' @param sealifebase_species_path Path to the SeaLifeBase `species.parquet`
#'   file.
#' @param sealifebase_ecology_path Path to the SeaLifeBase `ecology.parquet`
#'   file.
#' @param sealifebase_genera_path Path to the SeaLifeBase `genera.parquet`
#'   file.
#' @param sealifebase_families_path Path to the SeaLifeBase `families.parquet`
#'   file.
#' @param sealifebase_orders_path Path to the SeaLifeBase `orders.parquet`
#'   file.
#' @param sealifebase_classes_path Path to the SeaLifeBase `classes.parquet`
#'   file.
#' @param sealifebase_phylums_path Path to the SeaLifeBase `phylums.parquet`
#'   file.
#' @param sealifebase_phylumclass_path Path to the SeaLifeBase
#'   `phylumclass.parquet` file.
#' @param download_if_missing Logical scalar. If `TRUE`, missing default parquet
#'   files are downloaded before reading.
#' @param include_higher_ranks Logical scalar. If `TRUE`, also build
#'   standard-rank ecology evidence rows by aggregating the species-level table.
#' @param output_path Optional path to save the cleaned table as an RDS file.
#'
#' @return A list with `species`, `higher_ranks`, `combined`, and
#'   `diagnostics`. `species` contains one row per source/species `spec_code`;
#'   `higher_ranks` contains aggregated kingdom/phylum/class/order/family/genus
#'   rows when `include_higher_ranks = TRUE`.
#' @export
#' @seealso [build_fishbase_sealifebase_higher_rank_ecology_db()]
#'
#' @examples
#' \dontrun{
#' db <- build_fishbase_sealifebase_ecology_db(
#'   output_path = "data-raw/fishbase_sealifebase_ecology_db.rds"
#' )
#'
#' dplyr::count(db$species, source)
#' dplyr::count(db$higher_ranks, source, taxon_rank)
#' }
build_fishbase_sealifebase_ecology_db <- function(
  fishbase_species_path = "data-raw/fishbase_parquet/fb/species.parquet",
  fishbase_ecology_path = "data-raw/fishbase_parquet/fb/ecology.parquet",
  fishbase_genera_path = "data-raw/fishbase_parquet/fb/genera.parquet",
  fishbase_families_path = "data-raw/fishbase_parquet/fb/families.parquet",
  fishbase_orders_path = "data-raw/fishbase_parquet/fb/orders.parquet",
  fishbase_classes_path = "data-raw/fishbase_parquet/fb/classes.parquet",
  sealifebase_species_path = "data-raw/fishbase_parquet/slb/species.parquet",
  sealifebase_ecology_path = "data-raw/fishbase_parquet/slb/ecology.parquet",
  sealifebase_genera_path = "data-raw/fishbase_parquet/slb/genera.parquet",
  sealifebase_families_path = "data-raw/fishbase_parquet/slb/families.parquet",
  sealifebase_orders_path = "data-raw/fishbase_parquet/slb/orders.parquet",
  sealifebase_classes_path = "data-raw/fishbase_parquet/slb/classes.parquet",
  sealifebase_phylums_path = "data-raw/fishbase_parquet/slb/phylums.parquet",
  sealifebase_phylumclass_path = "data-raw/fishbase_parquet/slb/phylumclass.parquet",
  download_if_missing = TRUE,
  include_higher_ranks = TRUE,
  output_path = NULL
) {
  paths <- c(
    fishbase_species_path = .fb_required_scalar_character(fishbase_species_path, "fishbase_species_path"),
    fishbase_ecology_path = .fb_required_scalar_character(fishbase_ecology_path, "fishbase_ecology_path"),
    fishbase_genera_path = .fb_required_scalar_character(fishbase_genera_path, "fishbase_genera_path"),
    fishbase_families_path = .fb_required_scalar_character(fishbase_families_path, "fishbase_families_path"),
    fishbase_orders_path = .fb_required_scalar_character(fishbase_orders_path, "fishbase_orders_path"),
    fishbase_classes_path = .fb_required_scalar_character(fishbase_classes_path, "fishbase_classes_path"),
    sealifebase_species_path = .fb_required_scalar_character(sealifebase_species_path, "sealifebase_species_path"),
    sealifebase_ecology_path = .fb_required_scalar_character(sealifebase_ecology_path, "sealifebase_ecology_path"),
    sealifebase_genera_path = .fb_required_scalar_character(sealifebase_genera_path, "sealifebase_genera_path"),
    sealifebase_families_path = .fb_required_scalar_character(sealifebase_families_path, "sealifebase_families_path"),
    sealifebase_orders_path = .fb_required_scalar_character(sealifebase_orders_path, "sealifebase_orders_path"),
    sealifebase_classes_path = .fb_required_scalar_character(sealifebase_classes_path, "sealifebase_classes_path"),
    sealifebase_phylums_path = .fb_required_scalar_character(sealifebase_phylums_path, "sealifebase_phylums_path"),
    sealifebase_phylumclass_path = .fb_required_scalar_character(sealifebase_phylumclass_path, "sealifebase_phylumclass_path")
  )
  download_if_missing <- .fb_logical_scalar(download_if_missing, "download_if_missing")
  include_higher_ranks <- .fb_logical_scalar(include_higher_ranks, "include_higher_ranks")
  output_path <- .fb_optional_scalar_character(output_path, "output_path")

  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0 && isTRUE(download_if_missing)) {
    output_dir <- dirname(dirname(unname(paths[["fishbase_species_path"]])))
    download_fishbase_sealifebase_parquet(output_dir = output_dir)
    missing_paths <- paths[!file.exists(paths)]
  }
  if (length(missing_paths) > 0) {
    stop(
      "Missing FishBase/SeaLifeBase parquet file(s): ",
      paste(unname(missing_paths), collapse = ", "),
      ". Run `download_fishbase_sealifebase_parquet()` or download them ",
      "manually with command-line curl.",
      call. = FALSE
    )
  }

  .require_arrow_for_fishbase_parquet()
  out <- .combine_fishbase_sealifebase_ecology_tables(
    fishbase_species = arrow::read_parquet(paths[["fishbase_species_path"]]),
    fishbase_ecology = arrow::read_parquet(paths[["fishbase_ecology_path"]]),
    fishbase_genera = arrow::read_parquet(paths[["fishbase_genera_path"]]),
    fishbase_families = arrow::read_parquet(paths[["fishbase_families_path"]]),
    fishbase_orders = arrow::read_parquet(paths[["fishbase_orders_path"]]),
    fishbase_classes = arrow::read_parquet(paths[["fishbase_classes_path"]]),
    sealifebase_species = arrow::read_parquet(paths[["sealifebase_species_path"]]),
    sealifebase_ecology = arrow::read_parquet(paths[["sealifebase_ecology_path"]]),
    sealifebase_genera = arrow::read_parquet(paths[["sealifebase_genera_path"]]),
    sealifebase_families = arrow::read_parquet(paths[["sealifebase_families_path"]]),
    sealifebase_orders = arrow::read_parquet(paths[["sealifebase_orders_path"]]),
    sealifebase_classes = arrow::read_parquet(paths[["sealifebase_classes_path"]]),
    sealifebase_phylums = arrow::read_parquet(paths[["sealifebase_phylums_path"]]),
    sealifebase_phylumclass = arrow::read_parquet(paths[["sealifebase_phylumclass_path"]])
  )
  species <- out$species
  higher_ranks <- if (isTRUE(include_higher_ranks)) {
    build_fishbase_sealifebase_higher_rank_ecology_db(species)
  } else {
    tibble::as_tibble(data.frame())
  }
  combined <- .combine_species_and_higher_rank_fishbase_db(species, higher_ranks)
  result <- list(
    species = species,
    higher_ranks = higher_ranks,
    combined = combined,
    diagnostics = out$diagnostics
  )

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, output_path)
  }

  result
}

#' Build higher-rank FishBase/SeaLifeBase ecology evidence rows
#'
#' @description
#' Aggregates a species-level FishBase/SeaLifeBase ecology database to standard
#' kingdom, phylum, class, order, family, and genus rows. Intermediate ranks such
#' as superclass, subclass, and subfamily can remain as lineage metadata on
#' species rows, but they are not emitted as `taxon_rank` rows. Aggregated
#' ecology fields are derived from species-level evidence only. Native
#' rank-table metadata, when available on the `species_ecology_db` object as a
#' `"rank_metadata"` attribute, is retained in metadata columns and does not
#' override aggregated ecology fields.
#'
#' @param species_ecology_db Species-level table returned in the `species`
#'   element of [build_fishbase_sealifebase_ecology_db()].
#'
#' @return A tibble with aggregated higher-rank ecology evidence rows.
#' @export
#'
#' @examples
#' \dontrun{
#' db <- build_fishbase_sealifebase_ecology_db()
#' higher <- build_fishbase_sealifebase_higher_rank_ecology_db(db$species)
#' dplyr::count(higher, source, taxon_rank)
#' }
build_fishbase_sealifebase_higher_rank_ecology_db <- function(species_ecology_db) {
  if (!inherits(species_ecology_db, "data.frame")) {
    stop("`species_ecology_db` must be a data.frame.", call. = FALSE)
  }
  species_ecology_db <- tibble::as_tibble(species_ecology_db)
  required <- c("source", "taxon_rank", "taxon_name", "genus", "family", "order", "class", "phylum", "kingdom")
  missing <- setdiff(required, colnames(species_ecology_db))
  if (length(missing) > 0) {
    stop(
      "`species_ecology_db` is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  rank_rows <- lapply(.fishbase_standard_higher_ranks(), function(rank) {
    .aggregate_fishbase_rank(species_ecology_db, rank)
  })
  out <- dplyr::bind_rows(rank_rows)
  metadata <- attr(species_ecology_db, "rank_metadata", exact = TRUE)
  if (is.null(metadata)) {
    metadata <- .empty_fishbase_rank_metadata()
  }
  out <- dplyr::left_join(out, metadata,
    by = c("source", "taxon_rank", "taxon_name"),
    suffix = c("", "_metadata")
  )
  if ("common_name_metadata" %in% colnames(out)) {
    out$common_name <- .fb_coalesce_character(out$common_name_metadata, out$common_name)
    out$common_name_metadata <- NULL
  }
  .order_fishbase_higher_rank_columns(out)
}

.fishbase_standard_higher_ranks <- function() {
  c("kingdom", "phylum", "class", "order", "family", "genus")
}

.fishbase_standard_combined_ranks <- function() {
  c(.fishbase_standard_higher_ranks(), "species")
}

.aggregate_fishbase_rank <- function(species_db, rank) {
  rank_values <- as.character(species_db[[rank]])
  use <- .fb_non_empty(rank_values)
  if (!any(use)) {
    return(tibble::as_tibble(data.frame()))
  }
  rows <- species_db[use, , drop = FALSE]
  rows$taxon_rank <- rank
  rows$taxon_name <- rank_values[use]
  rows <- data.table::as.data.table(rows)
  group_cols <- c("source", "taxon_rank", "taxon_name")
  summaries <- list(
    unique(rows[, group_cols, with = FALSE]),
    rows[, .(
      spec_code = NA_character_,
      species = NA_character_,
      common_name = NA_character_,
      raw_demers_pelag = NA_character_,
      species_comments = NA_character_,
      n_species = .N,
      depth_min = .fb_min_numeric(depth_min),
      depth_max = .fb_max_numeric(depth_max),
      common_depth_min = .fb_min_numeric(common_depth_min),
      common_depth_max = .fb_max_numeric(common_depth_max)
    ), by = .(source, taxon_rank, taxon_name)]
  )
  for (column_name in c(
    "kingdom", "phylum", "superclass", "class", "subclass", "order",
    "family", "subfamily", "genus", "taxonomy_source"
  )) {
    summaries[[length(summaries) + 1]] <- .fishbase_collapse_by_group(
      rows = rows,
      field = column_name,
      output_column = column_name,
      split = FALSE
    )
  }
  for (column_name in c(
    "ecology_environment", "ecology_position", "ecology_habitat_broad",
    "ecology_depth_zone", "ecology_zone", "ecology_water_column_zone",
    "ecology_mobility", "ecology_size_class", "ecology_substrate",
    "ecology_special_habitat", "ecology_comments"
  )) {
    summaries[[length(summaries) + 1]] <- .fishbase_collapse_by_group(
      rows = rows,
      field = column_name,
      output_column = column_name,
      split = TRUE
    )
  }
  summaries[[length(summaries) + 1]] <- .fishbase_count_by_group(rows, "ecology_environment", "n_species_with_environment")
  summaries[[length(summaries) + 1]] <- .fishbase_count_by_group(rows, "ecology_habitat_broad", "n_species_with_habitat_broad")
  summaries[[length(summaries) + 1]] <- .fishbase_count_by_group(rows, "ecology_depth_zone", "n_species_with_depth")
  summaries[[length(summaries) + 1]] <- .fishbase_count_by_group(rows, "ecology_zone", "n_species_with_ecology_zone")
  summaries[[length(summaries) + 1]] <- .fishbase_count_by_group(rows, "ecology_water_column_zone", "n_species_with_water_column_zone")
  summaries[[length(summaries) + 1]] <- .fishbase_dominant_by_group(rows, "ecology_environment", "dominant_environment", "dominant_environment_prop")
  summaries[[length(summaries) + 1]] <- .fishbase_dominant_by_group(rows, "ecology_habitat_broad", "dominant_habitat_broad", "dominant_habitat_broad_prop")
  summaries[[length(summaries) + 1]] <- .fishbase_dominant_by_group(rows, "ecology_depth_zone", "dominant_depth_zone", "dominant_depth_zone_prop")

  out <- Reduce(function(left, right) {
    merge(left, right, by = group_cols, all.x = TRUE, sort = FALSE)
  }, summaries)
  for (column_name in c(
    "n_species_with_environment", "n_species_with_habitat_broad",
    "n_species_with_depth", "n_species_with_ecology_zone",
    "n_species_with_water_column_zone"
  )) {
    out[[column_name]][is.na(out[[column_name]])] <- 0L
  }
  out$ecology_evidence_strength <- vapply(seq_len(nrow(out)), function(index) {
    .fishbase_evidence_strength(
      n_species_with_habitat_broad = out$n_species_with_habitat_broad[[index]],
      dominant_habitat_broad_prop = out$dominant_habitat_broad_prop[[index]],
      ecology_habitat_broad = out$ecology_habitat_broad[[index]]
    )
  }, character(1))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out <- .fb_trim_lineage_below_rank_table(out, rank)
  tibble::as_tibble(out)
}

.fishbase_group_columns <- function() {
  c("source", "taxon_rank", "taxon_name")
}

.fishbase_tokens_by_group <- function(rows, field, split = TRUE, keep_unknown = FALSE) {
  group_cols <- .fishbase_group_columns()
  out <- rows[, group_cols, with = FALSE]
  out$.fishbase_row_id <- seq_len(nrow(rows))
  out$value <- trimws(as.character(rows[[field]]))
  out <- out[.fb_non_empty(value)]
  if (!isTRUE(keep_unknown)) {
    out <- out[tolower(value) != "unknown"]
  }
  if (nrow(out) == 0) {
    return(data.table::data.table(
      source = character(),
      taxon_rank = character(),
      taxon_name = character(),
      value = character()
    ))
  }
  if (isTRUE(split) && any(grepl(";", out$value, fixed = TRUE))) {
    out <- out[, .(
      value = trimws(unlist(strsplit(value, "\\s*;\\s*"), use.names = FALSE))
    ), by = .(source, taxon_rank, taxon_name, .fishbase_row_id)]
    out <- out[.fb_non_empty(value)]
    if (!isTRUE(keep_unknown)) {
      out <- out[tolower(value) != "unknown"]
    }
    out <- unique(out, by = c("source", "taxon_rank", "taxon_name", ".fishbase_row_id", "value"))
  }
  out[, .fishbase_row_id := NULL]
  out
}

.fishbase_collapse_by_group <- function(rows, field, output_column, split = TRUE) {
  tokens <- .fishbase_tokens_by_group(rows, field, split = split, keep_unknown = TRUE)
  if (nrow(tokens) == 0) {
    out <- data.table::data.table(
      source = character(),
      taxon_rank = character(),
      taxon_name = character(),
      value = character()
    )
  } else {
    out <- tokens[, .(
      value = paste(sort(unique(value)), collapse = "; ")
    ), by = .(source, taxon_rank, taxon_name)]
  }
  data.table::setnames(out, "value", output_column)
  out
}

.fishbase_count_by_group <- function(rows, field, output_column) {
  group_cols <- .fishbase_group_columns()
  out <- rows[, group_cols, with = FALSE]
  out$value <- trimws(as.character(rows[[field]]))
  out <- out[.fb_non_empty(value) & tolower(value) != "unknown"]
  if (nrow(out) == 0) {
    result <- data.table::data.table(
      source = character(),
      taxon_rank = character(),
      taxon_name = character(),
      n = integer()
    )
  } else {
    result <- out[, .(n = .N), by = .(source, taxon_rank, taxon_name)]
  }
  data.table::setnames(result, "n", output_column)
  result
}

.fishbase_dominant_by_group <- function(rows, field, value_column, prop_column) {
  tokens <- .fishbase_tokens_by_group(rows, field, split = TRUE, keep_unknown = TRUE)
  if (nrow(tokens) == 0) {
    out <- data.table::data.table(
      source = character(),
      taxon_rank = character(),
      taxon_name = character(),
      value = character(),
      prop = numeric()
    )
  } else {
    counts <- tokens[, .N, by = .(source, taxon_rank, taxon_name, value)]
    totals <- counts[, .(total = sum(N)), by = .(source, taxon_rank, taxon_name)]
    data.table::setorder(counts, source, taxon_rank, taxon_name, -N, value)
    out <- counts[, .SD[1], by = .(source, taxon_rank, taxon_name)]
    out <- merge(out, totals, by = c("source", "taxon_rank", "taxon_name"), all.x = TRUE, sort = FALSE)
    out$prop <- out$N / out$total
    out <- out[, .(source, taxon_rank, taxon_name, value, prop)]
  }
  data.table::setnames(out, c("value", "prop"), c(value_column, prop_column))
  out
}

.summarise_fishbase_rank_group <- function(group) {
  first_value <- function(column_name) {
    .collapse_fishbase_text_values(group[[column_name]])
  }
  environment_dom <- .dominant_fishbase_value(group$ecology_environment)
  habitat_dom <- .dominant_fishbase_value(group$ecology_habitat_broad)
  depth_dom <- .dominant_fishbase_value(group$ecology_depth_zone)

  data.frame(
    source = group$source[[1]],
    spec_code = NA_character_,
    taxon_rank = group$taxon_rank[[1]],
    taxon_name = group$taxon_name[[1]],
    kingdom = first_value("kingdom"),
    phylum = first_value("phylum"),
    superclass = first_value("superclass"),
    class = first_value("class"),
    subclass = first_value("subclass"),
    order = first_value("order"),
    family = first_value("family"),
    subfamily = first_value("subfamily"),
    genus = first_value("genus"),
    species = NA_character_,
    common_name = NA_character_,
    taxonomy_source = first_value("taxonomy_source"),
    ecology_environment = .collapse_fishbase_text_values(group$ecology_environment),
    ecology_position = .collapse_fishbase_text_values(group$ecology_position),
    ecology_habitat_broad = .collapse_fishbase_text_values(group$ecology_habitat_broad),
    ecology_depth_zone = .collapse_fishbase_text_values(group$ecology_depth_zone),
    depth_min = suppressWarnings(min(as.numeric(group$depth_min), na.rm = TRUE)),
    depth_max = suppressWarnings(max(as.numeric(group$depth_max), na.rm = TRUE)),
    common_depth_min = suppressWarnings(min(as.numeric(group$common_depth_min), na.rm = TRUE)),
    common_depth_max = suppressWarnings(max(as.numeric(group$common_depth_max), na.rm = TRUE)),
    ecology_zone = .collapse_fishbase_text_values(group$ecology_zone),
    ecology_water_column_zone = .collapse_fishbase_text_values(group$ecology_water_column_zone),
    ecology_mobility = .collapse_fishbase_text_values(group$ecology_mobility),
    ecology_size_class = .collapse_fishbase_text_values(group$ecology_size_class),
    ecology_substrate = .collapse_fishbase_text_values(group$ecology_substrate),
    ecology_special_habitat = .collapse_fishbase_text_values(group$ecology_special_habitat),
    raw_demers_pelag = NA_character_,
    species_comments = NA_character_,
    ecology_comments = .collapse_fishbase_text_values(group$ecology_comments),
    n_species = nrow(group),
    n_species_with_environment = .count_fishbase_useful_values(group$ecology_environment),
    n_species_with_habitat_broad = .count_fishbase_useful_values(group$ecology_habitat_broad),
    n_species_with_depth = .count_fishbase_useful_values(group$ecology_depth_zone),
    n_species_with_ecology_zone = .count_fishbase_useful_values(group$ecology_zone),
    n_species_with_water_column_zone = .count_fishbase_useful_values(group$ecology_water_column_zone),
    dominant_environment = environment_dom$value,
    dominant_environment_prop = environment_dom$prop,
    dominant_habitat_broad = habitat_dom$value,
    dominant_habitat_broad_prop = habitat_dom$prop,
    dominant_depth_zone = depth_dom$value,
    dominant_depth_zone_prop = depth_dom$prop,
    ecology_evidence_strength = .fishbase_evidence_strength(
      n_species_with_habitat_broad = .count_fishbase_useful_values(group$ecology_habitat_broad),
      dominant_habitat_broad_prop = habitat_dom$prop,
      ecology_habitat_broad = .collapse_fishbase_text_values(group$ecology_habitat_broad)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    .fb_trim_lineage_below_rank() |>
    .fb_fix_infinite_depths()
}

.combine_species_and_higher_rank_fishbase_db <- function(species, higher_ranks) {
  if ("taxon_rank" %in% colnames(species)) {
    species <- species[tolower(as.character(species$taxon_rank)) == "species", , drop = FALSE]
  }
  if ("taxon_rank" %in% colnames(higher_ranks)) {
    higher_ranks <- higher_ranks[
      tolower(as.character(higher_ranks$taxon_rank)) %in% .fishbase_standard_higher_ranks(),
      ,
      drop = FALSE
    ]
  }
  if (nrow(higher_ranks) == 0) {
    return(.order_fishbase_combined_db_columns(species))
  }
  out <- dplyr::bind_rows(
    .order_fishbase_combined_db_columns(species),
    .order_fishbase_combined_db_columns(higher_ranks)
  )
  if ("taxon_rank" %in% colnames(out)) {
    out <- out[
      tolower(as.character(out$taxon_rank)) %in% .fishbase_standard_combined_ranks(),
      ,
      drop = FALSE
    ]
  }
  tibble::as_tibble(out)
}

.order_fishbase_higher_rank_columns <- function(x) {
  expected <- c(
    .fishbase_species_db_columns(),
    "n_species",
    "n_species_with_environment",
    "n_species_with_habitat_broad",
    "n_species_with_depth",
    "n_species_with_ecology_zone",
    "n_species_with_water_column_zone",
    "dominant_environment",
    "dominant_environment_prop",
    "dominant_habitat_broad",
    "dominant_habitat_broad_prop",
    "dominant_depth_zone",
    "dominant_depth_zone_prop",
    "ecology_evidence_strength",
    "rank_water_salinity",
    "rank_body_shape",
    "rank_habitat_note",
    "rank_distribution",
    "rank_diagnosis",
    "rank_remarks",
    "rank_etymology"
  )
  for (column_name in expected) {
    if (!(column_name %in% colnames(x))) {
      x[[column_name]] <- NA_character_
    }
  }
  tibble::as_tibble(x[, expected, drop = FALSE])
}

.order_fishbase_combined_db_columns <- function(x) {
  expected <- unique(c(.fishbase_species_db_columns(), .fishbase_higher_rank_extra_columns()))
  for (column_name in expected) {
    if (!(column_name %in% colnames(x))) {
      x[[column_name]] <- NA
    }
  }
  tibble::as_tibble(x[, expected, drop = FALSE])
}

.fishbase_higher_rank_extra_columns <- function() {
  c(
    "n_species",
    "n_species_with_environment",
    "n_species_with_habitat_broad",
    "n_species_with_depth",
    "n_species_with_ecology_zone",
    "n_species_with_water_column_zone",
    "dominant_environment",
    "dominant_environment_prop",
    "dominant_habitat_broad",
    "dominant_habitat_broad_prop",
    "dominant_depth_zone",
    "dominant_depth_zone_prop",
    "ecology_evidence_strength",
    "rank_water_salinity",
    "rank_body_shape",
    "rank_habitat_note",
    "rank_distribution",
    "rank_diagnosis",
    "rank_remarks",
    "rank_etymology"
  )
}

.fb_trim_lineage_below_rank <- function(x) {
  rank <- x$taxon_rank[[1]]
  keep <- switch(rank,
    genus = c("kingdom", "phylum", "class", "order", "family", "genus"),
    family = c("kingdom", "phylum", "class", "order", "family"),
    order = c("kingdom", "phylum", "class", "order"),
    class = c("kingdom", "phylum", "class"),
    phylum = c("kingdom", "phylum"),
    kingdom = c("kingdom"),
    c("kingdom", "phylum", "class", "order", "family", "genus")
  )
  lineage_columns <- c("kingdom", "phylum", "superclass", "class", "subclass", "order", "family", "subfamily", "genus", "species")
  drop_columns <- setdiff(lineage_columns, keep)
  for (column_name in drop_columns) {
    if (column_name %in% colnames(x)) {
      x[[column_name]] <- NA_character_
    }
  }
  x
}

.fb_trim_lineage_below_rank_table <- function(x, rank) {
  keep <- switch(rank,
    genus = c("kingdom", "phylum", "class", "order", "family", "genus"),
    family = c("kingdom", "phylum", "class", "order", "family"),
    order = c("kingdom", "phylum", "class", "order"),
    class = c("kingdom", "phylum", "class"),
    phylum = c("kingdom", "phylum"),
    kingdom = c("kingdom"),
    c("kingdom", "phylum", "class", "order", "family", "genus")
  )
  lineage_columns <- c("kingdom", "phylum", "superclass", "class", "subclass", "order", "family", "subfamily", "genus", "species")
  drop_columns <- setdiff(lineage_columns, keep)
  for (column_name in drop_columns) {
    if (column_name %in% colnames(x)) {
      x[[column_name]] <- NA_character_
    }
  }
  x
}

.fb_fix_infinite_depths <- function(x) {
  for (column_name in c("depth_min", "depth_max", "common_depth_min", "common_depth_max")) {
    if (column_name %in% colnames(x)) {
      value <- suppressWarnings(as.numeric(x[[column_name]]))
      value[is.infinite(value)] <- NA_real_
      x[[column_name]] <- value
    }
  }
  x
}

.fb_min_numeric <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (all(is.na(value))) {
    return(NA_real_)
  }
  min(value, na.rm = TRUE)
}

.fb_max_numeric <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (all(is.na(value))) {
    return(NA_real_)
  }
  max(value, na.rm = TRUE)
}

.count_fishbase_useful_values <- function(values) {
  values <- trimws(as.character(values))
  values <- values[.fb_non_empty(values)]
  sum(tolower(values) != "unknown")
}

.dominant_fishbase_value <- function(values) {
  values <- trimws(as.character(values))
  values <- values[.fb_non_empty(values)]
  if (length(values) == 0) {
    return(list(value = NA_character_, prop = NA_real_))
  }
  split <- unlist(lapply(values, function(value) {
    terms <- unlist(strsplit(value, "\\s*;\\s*"), use.names = FALSE)
    terms <- trimws(terms)
    unique(terms[.fb_non_empty(terms)])
  }), use.names = FALSE)
  if (length(split) == 0) {
    return(list(value = NA_character_, prop = NA_real_))
  }
  counts <- as.data.frame(table(split), stringsAsFactors = FALSE)
  colnames(counts) <- c("value", "n")
  counts <- counts[order(-counts$n, counts$value), , drop = FALSE]
  list(
    value = as.character(counts$value[[1]]),
    prop = as.numeric(counts$n[[1]]) / sum(counts$n)
  )
}

.fb_useful_split_values <- function(value) {
  split <- .split_fishbase_collapsed_text(value)
  split <- split[.fb_non_empty(split)]
  split <- split[tolower(split) != "unknown"]
  split
}

.fishbase_evidence_strength <- function(n_species_with_habitat_broad, dominant_habitat_broad_prop, ecology_habitat_broad) {
  habitats <- .fb_useful_split_values(ecology_habitat_broad)
  if (length(habitats) == 0 || n_species_with_habitat_broad == 0) {
    return("unknown")
  }
  if (length(unique(habitats)) > 1 && (is.na(dominant_habitat_broad_prop) || dominant_habitat_broad_prop < 0.6)) {
    return("mixed")
  }
  if (n_species_with_habitat_broad >= 10 && !is.na(dominant_habitat_broad_prop) && dominant_habitat_broad_prop >= 0.8) {
    return("high")
  }
  if (n_species_with_habitat_broad >= 5 && !is.na(dominant_habitat_broad_prop) && dominant_habitat_broad_prop >= 0.6) {
    return("moderate")
  }
  "low"
}

.fishbase_sealifebase_parquet_manifest <- function(output_dir, version) {
  base_url <- paste0("https://data.source.coop/cboettig/fishbase")
  rows <- data.frame(
    source = c(rep("fishbase", 6), rep("sealifebase", 8)),
    source_dir = c(rep("fb", 6), rep("slb", 8)),
    table = c(
      "species", "ecology", "genera", "families", "orders", "classes",
      "species", "ecology", "genera", "families", "orders", "classes",
      "phylums", "phylumclass"
    ),
    stringsAsFactors = FALSE
  )
  rows$url <- paste0(
    base_url,
    "/",
    rows$source_dir,
    "/v",
    version,
    "/parquet/",
    rows$table,
    ".parquet"
  )
  rows$path <- file.path(output_dir, rows$source_dir, paste0(rows$table, ".parquet"))
  rows$downloaded <- NA
  rows[, c("source", "table", "url", "path", "downloaded"), drop = FALSE]
}

.combine_fishbase_sealifebase_ecology_tables <- function(
  fishbase_species,
  fishbase_ecology,
  fishbase_genera = NULL,
  fishbase_families = NULL,
  fishbase_orders = NULL,
  fishbase_classes = NULL,
  sealifebase_species,
  sealifebase_ecology,
  sealifebase_genera = NULL,
  sealifebase_families = NULL,
  sealifebase_orders = NULL,
  sealifebase_classes = NULL,
  sealifebase_phylums = NULL,
  sealifebase_phylumclass = NULL
) {
  fishbase <- .join_fishbase_species_ecology(
    species = fishbase_species,
    ecology = fishbase_ecology,
    taxonomy = clean_fishbase_taxonomy_table(
      species = fishbase_species,
      genera = fishbase_genera,
      families = fishbase_families,
      orders = fishbase_orders,
      classes = fishbase_classes
    ),
    source = "fishbase"
  )
  fishbase_report <- attr(fishbase, "build_report", exact = TRUE)
  sealifebase <- .join_fishbase_species_ecology(
    species = sealifebase_species,
    ecology = sealifebase_ecology,
    taxonomy = clean_sealifebase_taxonomy_table(
      species = sealifebase_species,
      genera = sealifebase_genera,
      families = sealifebase_families,
      orders = sealifebase_orders,
      classes = sealifebase_classes,
      phylums = sealifebase_phylums,
      phylumclass = sealifebase_phylumclass
    ),
    source = "sealifebase"
  )
  sealifebase_report <- attr(sealifebase, "build_report", exact = TRUE)

  out <- .order_fishbase_ecology_db_columns(tibble::as_tibble(
    dplyr::bind_rows(fishbase, sealifebase)
  ))
  expected_rows <- nrow(as.data.frame(fishbase_species)) +
    nrow(as.data.frame(sealifebase_species))
  actual_rows <- nrow(out)
  report <- data.frame(
    fishbase_species_rows = fishbase_report$species_rows,
    fishbase_ecology_rows = fishbase_report$ecology_rows,
    fishbase_ecology_distinct_spec_code = fishbase_report$ecology_distinct_spec_code,
    fishbase_ecology_unmatched_spec_code = fishbase_report$ecology_unmatched_spec_code,
    fishbase_taxonomy_rows = fishbase_report$taxonomy_rows,
    fishbase_taxonomy_unmatched_spec_code = fishbase_report$taxonomy_unmatched_spec_code,
    sealifebase_species_rows = sealifebase_report$species_rows,
    sealifebase_ecology_rows = sealifebase_report$ecology_rows,
    sealifebase_ecology_distinct_spec_code = sealifebase_report$ecology_distinct_spec_code,
    sealifebase_ecology_unmatched_spec_code = sealifebase_report$ecology_unmatched_spec_code,
    sealifebase_taxonomy_rows = sealifebase_report$taxonomy_rows,
    sealifebase_taxonomy_unmatched_spec_code = sealifebase_report$taxonomy_unmatched_spec_code,
    expected_rows = expected_rows,
    final_rows = actual_rows,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (actual_rows != expected_rows) {
    stop(
      "FishBase/SeaLifeBase ecology DB row-count check failed: expected ",
      expected_rows,
      " species rows but built ",
      actual_rows,
      " rows. This usually means ecology rows were not collapsed before joining.",
      call. = FALSE
    )
  }
  attr(out, "build_report") <- report
  attr(out, "rank_metadata") <- dplyr::bind_rows(
    attr(fishbase, "rank_metadata", exact = TRUE),
    attr(sealifebase, "rank_metadata", exact = TRUE)
  )
  list(species = out, diagnostics = report)
}

.join_fishbase_species_ecology <- function(species, ecology, taxonomy, source) {
  species_rows <- nrow(as.data.frame(species))
  ecology_rows <- nrow(as.data.frame(ecology))
  clean_species <- clean_fishbase_species_table(species, source = source)
  clean_ecology <- clean_fishbase_ecology_table(ecology, source = source)
  clean_taxonomy <- .fb_taxonomy_table_or_empty(taxonomy, source = source)
  out <- dplyr::left_join(clean_species, clean_ecology, by = c("source", "spec_code"))
  out <- dplyr::left_join(out, clean_taxonomy, by = c("source", "spec_code"), suffix = c("", "_taxonomy"))
  out <- .coalesce_joined_taxonomy_columns(out)
  out <- .order_fishbase_ecology_db_columns(out)
  expected_rows <- nrow(clean_species)
  if (nrow(out) != expected_rows) {
    stop(
      "FishBase/SeaLifeBase ecology DB row-count check failed for source `",
      source,
      "`: expected ",
      expected_rows,
      " species rows but built ",
      nrow(out),
      " rows.",
      call. = FALSE
    )
  }
  attr(out, "build_report") <- data.frame(
    source = source,
    species_rows = species_rows,
    ecology_rows = ecology_rows,
    ecology_distinct_spec_code = length(unique(clean_ecology$spec_code[.fb_non_empty(clean_ecology$spec_code)])),
    ecology_unmatched_spec_code = sum(
      !(clean_ecology$spec_code[.fb_non_empty(clean_ecology$spec_code)] %in% clean_species$spec_code)
    ),
    taxonomy_rows = nrow(clean_taxonomy),
    taxonomy_unmatched_spec_code = sum(
      !(clean_taxonomy$spec_code[.fb_non_empty(clean_taxonomy$spec_code)] %in% clean_species$spec_code)
    ),
    final_rows = nrow(out),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  attr(out, "rank_metadata") <- attr(clean_taxonomy, "rank_metadata", exact = TRUE)
  out
}

clean_fishbase_species_table <- function(species, source) {
  source <- .fb_required_scalar_character(source, "source")
  species <- .clean_fishbase_names(species, "species")
  species <- .fb_add_missing_columns(
    species,
    c(
      "spec_code",
      "genus",
      "species",
      "fbname",
      "f_bname",
      "fresh",
      "brack",
      "saltwater",
      "land",
      "demers_pelag",
      "depth_range_shallow",
      "depth_range_deep",
      "depth_range_com_shallow",
      "depth_range_com_deep",
      "comments"
    )
  )

  genus <- as.character(species$genus)
  species_name <- as.character(species$species)
  common_name <- .fb_first_non_empty_column(species, c("fbname", "f_bname", "common_name"))
  taxon_name <- ifelse(
    .fb_non_empty(genus) & .fb_non_empty(species_name),
    paste(genus, species_name),
    NA_character_
  )
  ecology_position <- .fishbase_position_text(species$demers_pelag)

  out <- data.frame(
    source = source,
    spec_code = .fb_as_spec_code(species$spec_code),
    taxon_rank = "species",
    taxon_name = taxon_name,
    genus = genus,
    species = species_name,
    common_name = common_name,
    ecology_environment = .fishbase_environment_text(species),
    ecology_position = ecology_position,
    ecology_habitat_broad = ecology_position,
    ecology_depth_zone = .fishbase_depth_zone(species$depth_range_deep),
    depth_min = suppressWarnings(as.numeric(species$depth_range_shallow)),
    depth_max = suppressWarnings(as.numeric(species$depth_range_deep)),
    common_depth_min = suppressWarnings(as.numeric(species$depth_range_com_shallow)),
    common_depth_max = suppressWarnings(as.numeric(species$depth_range_com_deep)),
    raw_demers_pelag = as.character(species$demers_pelag),
    species_comments = as.character(species$comments),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  tibble::as_tibble(out)
}

clean_fishbase_ecology_table <- function(ecology, source) {
  source <- .fb_required_scalar_character(source, "source")
  ecology <- .clean_fishbase_names(ecology, "ecology")
  ecology <- .fb_add_missing_columns(
    ecology,
    unique(c(
      "spec_code",
      "add_rems",
      names(.fishbase_ecology_zone_map()),
      names(.fishbase_water_column_zone_map()),
      names(.fishbase_mobility_map()),
      names(.fishbase_size_class_map()),
      names(.fishbase_substrate_map()),
      names(.fishbase_special_habitat_map())
    ))
  )

  out <- data.frame(
    source = source,
    spec_code = .fb_as_spec_code(ecology$spec_code),
    ecology_comments = as.character(ecology$add_rems),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[.fb_non_empty(out$spec_code), , drop = FALSE]
  if (nrow(out) == 0) {
    return(.empty_clean_fishbase_ecology_table())
  }

  out <- unique(out)
  out <- .collapse_fishbase_ecology_comments(out)
  collapsed_tables <- list(
    .collapse_fishbase_flag_table(ecology, source, .fishbase_ecology_zone_map(), "ecology_zone"),
    .collapse_fishbase_flag_table(ecology, source, .fishbase_water_column_zone_map(), "ecology_water_column_zone"),
    .collapse_fishbase_flag_table(ecology, source, .fishbase_mobility_map(), "ecology_mobility"),
    .collapse_fishbase_flag_table(ecology, source, .fishbase_size_class_map(), "ecology_size_class"),
    .collapse_fishbase_flag_table(ecology, source, .fishbase_substrate_map(), "ecology_substrate"),
    .collapse_fishbase_flag_table(ecology, source, .fishbase_special_habitat_map(), "ecology_special_habitat")
  )
  for (collapsed_table in collapsed_tables) {
    out <- merge(out, collapsed_table,
      by = c("source", "spec_code"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  .order_clean_fishbase_ecology_columns(tibble::as_tibble(out))
}

.collapse_fishbase_ecology_comments <- function(x) {
  comments <- x[, c("source", "spec_code", "ecology_comments"), drop = FALSE]
  comments <- comments[.fb_non_empty(comments$spec_code), , drop = FALSE]
  if (nrow(comments) == 0) {
    x$ecology_comments <- NA_character_
    return(unique(x[, c("source", "spec_code", "ecology_comments"), drop = FALSE]))
  }
  data.table::setDT(comments)
  comments <- comments[, .(
    ecology_comments = .collapse_fishbase_text_values(ecology_comments)
  ), by = .(source, spec_code)]
  as.data.frame(comments, stringsAsFactors = FALSE)
}

.empty_clean_fishbase_ecology_table <- function() {
  tibble::as_tibble(data.frame(
    source = character(),
    spec_code = character(),
    ecology_comments = character(),
    ecology_zone = character(),
    ecology_water_column_zone = character(),
    ecology_mobility = character(),
    ecology_size_class = character(),
    ecology_substrate = character(),
    ecology_special_habitat = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}

clean_fishbase_taxonomy_table <- function(
  species,
  genera = NULL,
  families = NULL,
  orders = NULL,
  classes = NULL
) {
  species <- .clean_fishbase_names(species, "fishbase_species_taxonomy")
  genera <- .clean_fishbase_names(.fb_null_data_frame(genera), "fishbase_genera")
  families <- .clean_fishbase_names(.fb_null_data_frame(families), "fishbase_families")
  orders <- .clean_fishbase_names(.fb_null_data_frame(orders), "fishbase_orders")
  classes <- .clean_fishbase_names(.fb_null_data_frame(classes), "fishbase_classes")

  species <- .fb_add_missing_columns(species, c("spec_code", "species", "genus", "subfamily", "gen_code", "fam_code"))
  metadata_columns <- c(
    "body_shape_i", "habitat", "water_salinity", "distribution",
    "diagnosis", "remark", "remarks", "comment", "classification_remark",
    "etymology", "etymology_1", "etymology1", "etymology_2", "etymology2",
    "name_etymology", "etym"
  )
  genera <- .fb_add_missing_columns(genera, c("gen_code", "gen_name", "gen_com_name", "fam_code", "subfamily", metadata_columns))
  families <- .fb_add_missing_columns(families, c("fam_code", "family", "common_name", "order", "ordnum", "class", "class_num", metadata_columns))
  orders <- .fb_add_missing_columns(orders, c("ordnum", "order", "common_name", "class_num", "class", metadata_columns))
  classes <- .fb_add_missing_columns(classes, c("class_num", "class", "common_name", "super_class", "subclass", metadata_columns))

  sp <- data.frame(
    source = "fishbase",
    spec_code = .fb_as_spec_code(species$spec_code),
    species_epithet = as.character(species$species),
    species_genus = as.character(species$genus),
    species_subfamily = as.character(species$subfamily),
    gen_code = .fb_as_spec_code(species$gen_code),
    fam_code = .fb_as_spec_code(species$fam_code),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  gen <- data.frame(
    gen_code = .fb_as_spec_code(genera$gen_code),
    genus_from_genera = as.character(genera$gen_name),
    genus_common_name = as.character(genera$gen_com_name),
    fam_code_from_genera = .fb_as_spec_code(genera$fam_code),
    subfamily_from_genera = as.character(genera$subfamily),
    body_shape_i = as.character(genera$body_shape_i),
    habitat = as.character(genera$habitat),
    water_salinity = as.character(genera$water_salinity),
    distribution = as.character(genera$distribution),
    diagnosis = as.character(genera$diagnosis),
    remark = as.character(genera$remark),
    remarks = as.character(genera$remarks),
    comment = as.character(genera$comment),
    classification_remark = as.character(genera$classification_remark),
    etymology = as.character(genera$etymology),
    etymology_1 = as.character(genera$etymology_1),
    etymology1 = as.character(genera$etymology1),
    etymology_2 = as.character(genera$etymology_2),
    etymology2 = as.character(genera$etymology2),
    name_etymology = as.character(genera$name_etymology),
    etym = as.character(genera$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  gen <- .fb_distinct_by_key(gen, "gen_code")
  fam <- data.frame(
    fam_code = .fb_as_spec_code(families$fam_code),
    family = as.character(families$family),
    family_common_name = as.character(families$common_name),
    order_from_family = as.character(families$order),
    ordnum = .fb_as_spec_code(families$ordnum),
    class_from_family = as.character(families$class),
    class_num = .fb_as_spec_code(families$class_num),
    body_shape_i = as.character(families$body_shape_i),
    habitat = as.character(families$habitat),
    water_salinity = as.character(families$water_salinity),
    distribution = as.character(families$distribution),
    diagnosis = as.character(families$diagnosis),
    remark = as.character(families$remark),
    remarks = as.character(families$remarks),
    comment = as.character(families$comment),
    classification_remark = as.character(families$classification_remark),
    etymology = as.character(families$etymology),
    etymology_1 = as.character(families$etymology_1),
    etymology1 = as.character(families$etymology1),
    etymology_2 = as.character(families$etymology_2),
    etymology2 = as.character(families$etymology2),
    name_etymology = as.character(families$name_etymology),
    etym = as.character(families$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fam <- .fb_distinct_by_key(fam, "fam_code")
  ord <- data.frame(
    ordnum = .fb_as_spec_code(orders$ordnum),
    order_from_orders = as.character(orders$order),
    order_common_name = as.character(orders$common_name),
    class_num_from_orders = .fb_as_spec_code(orders$class_num),
    class_from_orders = as.character(orders$class),
    body_shape_i = as.character(orders$body_shape_i),
    habitat = as.character(orders$habitat),
    water_salinity = as.character(orders$water_salinity),
    distribution = as.character(orders$distribution),
    diagnosis = as.character(orders$diagnosis),
    remark = as.character(orders$remark),
    remarks = as.character(orders$remarks),
    comment = as.character(orders$comment),
    classification_remark = as.character(orders$classification_remark),
    etymology = as.character(orders$etymology),
    etymology_1 = as.character(orders$etymology_1),
    etymology1 = as.character(orders$etymology1),
    etymology_2 = as.character(orders$etymology_2),
    etymology2 = as.character(orders$etymology2),
    name_etymology = as.character(orders$name_etymology),
    etym = as.character(orders$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ord <- .fb_distinct_by_key(ord, "ordnum")
  cls <- data.frame(
    class_num = .fb_as_spec_code(classes$class_num),
    class_from_classes = as.character(classes$class),
    class_common_name = as.character(classes$common_name),
    superclass = as.character(classes$super_class),
    subclass = as.character(classes$subclass),
    body_shape_i = as.character(classes$body_shape_i),
    habitat = as.character(classes$habitat),
    water_salinity = as.character(classes$water_salinity),
    distribution = as.character(classes$distribution),
    diagnosis = as.character(classes$diagnosis),
    remark = as.character(classes$remark),
    remarks = as.character(classes$remarks),
    comment = as.character(classes$comment),
    classification_remark = as.character(classes$classification_remark),
    etymology = as.character(classes$etymology),
    etymology_1 = as.character(classes$etymology_1),
    etymology1 = as.character(classes$etymology1),
    etymology_2 = as.character(classes$etymology_2),
    etymology2 = as.character(classes$etymology2),
    name_etymology = as.character(classes$name_etymology),
    etym = as.character(classes$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cls <- .fb_distinct_by_key(cls, "class_num")

  out <- dplyr::left_join(sp, gen, by = "gen_code")
  out$fam_code <- .fb_coalesce_character(out$fam_code, out$fam_code_from_genera)
  out <- dplyr::left_join(out, fam, by = "fam_code")
  out <- dplyr::left_join(out, ord, by = "ordnum")
  out$class_num <- .fb_coalesce_character(out$class_num, out$class_num_from_orders)
  out <- dplyr::left_join(out, cls, by = "class_num")
  out$genus <- .fb_coalesce_character(out$species_genus, out$genus_from_genera)
  out$species <- as.character(out$species_epithet)
  out$subfamily <- .fb_coalesce_character(out$species_subfamily, out$subfamily_from_genera)
  out$order <- .fb_coalesce_character(out$order_from_family, out$order_from_orders)
  out$class <- .fb_coalesce_character(out$class_from_classes, out$class_from_family, out$class_from_orders)
  out$kingdom <- "Metazoa"
  out$phylum <- "Chordata"
  out$taxonomy_source <- "fishbase_native_inferred_phylum"
  out$taxon_name <- ifelse(.fb_non_empty(out$genus) & .fb_non_empty(out$species), paste(out$genus, out$species), NA_character_)

  result <- .order_fishbase_taxonomy_columns(out)
  attr(result, "rank_metadata") <- dplyr::bind_rows(
    .rank_metadata_from_table(gen, "fishbase", "genus", "genus_from_genera", common_name_col = "genus_common_name"),
    .rank_metadata_from_table(fam, "fishbase", "family", "family", common_name_col = "family_common_name"),
    .rank_metadata_from_table(ord, "fishbase", "order", "order_from_orders", common_name_col = "order_common_name"),
    .rank_metadata_from_table(cls, "fishbase", "class", "class_from_classes", common_name_col = "class_common_name")
  )
  result
}

clean_sealifebase_taxonomy_table <- function(
  species,
  genera = NULL,
  families = NULL,
  orders = NULL,
  classes = NULL,
  phylums = NULL,
  phylumclass = NULL
) {
  species <- .clean_fishbase_names(species, "sealifebase_species_taxonomy")
  genera <- .clean_fishbase_names(.fb_null_data_frame(genera), "sealifebase_genera")
  families <- .clean_fishbase_names(.fb_null_data_frame(families), "sealifebase_families")
  orders <- .clean_fishbase_names(.fb_null_data_frame(orders), "sealifebase_orders")
  classes <- .clean_fishbase_names(.fb_null_data_frame(classes), "sealifebase_classes")
  phylums <- .clean_fishbase_names(.fb_null_data_frame(phylums), "sealifebase_phylums")
  phylumclass <- .clean_fishbase_names(.fb_null_data_frame(phylumclass), "sealifebase_phylumclass")

  species <- .fb_add_missing_columns(species, c("spec_code", "species", "genus", "gen_code", "fam_code"))
  metadata_columns <- c(
    "body_shape_i", "habitat", "water_salinity", "distribution",
    "diagnosis", "remark", "remarks", "comment", "classification_remark",
    "etymology", "etymology_1", "etymology1", "etymology_2", "etymology2",
    "name_etymology", "etym"
  )
  genera <- .fb_add_missing_columns(genera, c("gen_code", "gen_name", "common_name", "famcode", "fam_code", "subfamily", metadata_columns))
  families <- .fb_add_missing_columns(families, c("fam_code", "family", "common_name", "order", "ordnum", "class", "class_num", "phylum", metadata_columns))
  orders <- .fb_add_missing_columns(orders, c("ordnum", "order", "common_name", "class_num", "class", "phylum", "phylum_num", metadata_columns))
  classes <- .fb_add_missing_columns(classes, c("class_num", "class", "common_name", "phylum", "phylum_num", metadata_columns))
  phylums <- .fb_add_missing_columns(phylums, c("phylum_id", "kingdom", "phylum", "common_name", metadata_columns))
  phylumclass <- .fb_add_missing_columns(phylumclass, c("phylum", "class"))

  sp <- data.frame(
    source = "sealifebase",
    spec_code = .fb_as_spec_code(species$spec_code),
    species_epithet = as.character(species$species),
    species_genus = as.character(species$genus),
    gen_code = .fb_as_spec_code(species$gen_code),
    fam_code = .fb_as_spec_code(species$fam_code),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  gen <- data.frame(
    gen_code = .fb_as_spec_code(genera$gen_code),
    genus_from_genera = as.character(genera$gen_name),
    genus_common_name = as.character(genera$common_name),
    fam_code_from_genera = .fb_as_spec_code(.fb_coalesce_character(genera$fam_code, genera$famcode)),
    subfamily = as.character(genera$subfamily),
    body_shape_i = as.character(genera$body_shape_i),
    habitat = as.character(genera$habitat),
    water_salinity = as.character(genera$water_salinity),
    distribution = as.character(genera$distribution),
    diagnosis = as.character(genera$diagnosis),
    remark = as.character(genera$remark),
    remarks = as.character(genera$remarks),
    comment = as.character(genera$comment),
    classification_remark = as.character(genera$classification_remark),
    etymology = as.character(genera$etymology),
    etymology_1 = as.character(genera$etymology_1),
    etymology1 = as.character(genera$etymology1),
    etymology_2 = as.character(genera$etymology_2),
    etymology2 = as.character(genera$etymology2),
    name_etymology = as.character(genera$name_etymology),
    etym = as.character(genera$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  gen <- .fb_distinct_by_key(gen, "gen_code")
  fam <- data.frame(
    fam_code = .fb_as_spec_code(families$fam_code),
    family = as.character(families$family),
    family_common_name = as.character(families$common_name),
    order_from_family = as.character(families$order),
    ordnum = .fb_as_spec_code(families$ordnum),
    class_from_family = as.character(families$class),
    class_num = .fb_as_spec_code(families$class_num),
    phylum_from_family = as.character(families$phylum),
    body_shape_i = as.character(families$body_shape_i),
    habitat = as.character(families$habitat),
    water_salinity = as.character(families$water_salinity),
    distribution = as.character(families$distribution),
    diagnosis = as.character(families$diagnosis),
    remark = as.character(families$remark),
    remarks = as.character(families$remarks),
    comment = as.character(families$comment),
    classification_remark = as.character(families$classification_remark),
    etymology = as.character(families$etymology),
    etymology_1 = as.character(families$etymology_1),
    etymology1 = as.character(families$etymology1),
    etymology_2 = as.character(families$etymology_2),
    etymology2 = as.character(families$etymology2),
    name_etymology = as.character(families$name_etymology),
    etym = as.character(families$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fam <- .fb_distinct_by_key(fam, "fam_code")
  ord <- data.frame(
    ordnum = .fb_as_spec_code(orders$ordnum),
    order_from_orders = as.character(orders$order),
    order_common_name = as.character(orders$common_name),
    class_num_from_orders = .fb_as_spec_code(orders$class_num),
    class_from_orders = as.character(orders$class),
    phylum_from_orders = as.character(orders$phylum),
    phylum_id_from_orders = .fb_as_spec_code(orders$phylum_num),
    body_shape_i = as.character(orders$body_shape_i),
    habitat = as.character(orders$habitat),
    water_salinity = as.character(orders$water_salinity),
    distribution = as.character(orders$distribution),
    diagnosis = as.character(orders$diagnosis),
    remark = as.character(orders$remark),
    remarks = as.character(orders$remarks),
    comment = as.character(orders$comment),
    classification_remark = as.character(orders$classification_remark),
    etymology = as.character(orders$etymology),
    etymology_1 = as.character(orders$etymology_1),
    etymology1 = as.character(orders$etymology1),
    etymology_2 = as.character(orders$etymology_2),
    etymology2 = as.character(orders$etymology2),
    name_etymology = as.character(orders$name_etymology),
    etym = as.character(orders$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ord <- .fb_distinct_by_key(ord, "ordnum")
  cls <- data.frame(
    class_num = .fb_as_spec_code(classes$class_num),
    class_from_classes = as.character(classes$class),
    class_common_name = as.character(classes$common_name),
    phylum_from_classes = as.character(classes$phylum),
    phylum_id_from_classes = .fb_as_spec_code(classes$phylum_num),
    body_shape_i = as.character(classes$body_shape_i),
    habitat = as.character(classes$habitat),
    water_salinity = as.character(classes$water_salinity),
    distribution = as.character(classes$distribution),
    diagnosis = as.character(classes$diagnosis),
    remark = as.character(classes$remark),
    remarks = as.character(classes$remarks),
    comment = as.character(classes$comment),
    classification_remark = as.character(classes$classification_remark),
    etymology = as.character(classes$etymology),
    etymology_1 = as.character(classes$etymology_1),
    etymology1 = as.character(classes$etymology1),
    etymology_2 = as.character(classes$etymology_2),
    etymology2 = as.character(classes$etymology2),
    name_etymology = as.character(classes$name_etymology),
    etym = as.character(classes$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cls <- .fb_distinct_by_key(cls, "class_num")
  pc <- data.frame(
    phylum_from_phylumclass = as.character(phylumclass$phylum),
    class_from_phylumclass = as.character(phylumclass$class),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pc <- .fb_distinct_by_key(pc, "class_from_phylumclass")
  phy <- data.frame(
    phylum_id = .fb_as_spec_code(phylums$phylum_id),
    kingdom_from_phylums = as.character(phylums$kingdom),
    phylum_from_phylums = as.character(phylums$phylum),
    phylum_common_name = as.character(phylums$common_name),
    body_shape_i = as.character(phylums$body_shape_i),
    habitat = as.character(phylums$habitat),
    water_salinity = as.character(phylums$water_salinity),
    distribution = as.character(phylums$distribution),
    diagnosis = as.character(phylums$diagnosis),
    remark = as.character(phylums$remark),
    remarks = as.character(phylums$remarks),
    comment = as.character(phylums$comment),
    classification_remark = as.character(phylums$classification_remark),
    etymology = as.character(phylums$etymology),
    etymology_1 = as.character(phylums$etymology_1),
    etymology1 = as.character(phylums$etymology1),
    etymology_2 = as.character(phylums$etymology_2),
    etymology2 = as.character(phylums$etymology2),
    name_etymology = as.character(phylums$name_etymology),
    etym = as.character(phylums$etym),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  phy <- .fb_distinct_by_key(phy, "phylum_from_phylums")

  out <- dplyr::left_join(sp, gen, by = "gen_code")
  out$fam_code <- .fb_coalesce_character(out$fam_code, out$fam_code_from_genera)
  out <- dplyr::left_join(out, fam, by = "fam_code")
  out <- dplyr::left_join(out, ord, by = "ordnum")
  out$class_num <- .fb_coalesce_character(out$class_num, out$class_num_from_orders)
  out <- dplyr::left_join(out, cls, by = "class_num")
  out$class <- .fb_coalesce_character(out$class_from_classes, out$class_from_family, out$class_from_orders)
  out <- dplyr::left_join(out, pc, by = c("class" = "class_from_phylumclass"))
  out$phylum <- .fb_coalesce_character(
    out$phylum_from_family,
    out$phylum_from_orders,
    out$phylum_from_classes,
    out$phylum_from_phylumclass
  )
  out$phylum_id <- .fb_coalesce_character(out$phylum_id_from_classes, out$phylum_id_from_orders)
  out <- dplyr::left_join(out, phy, by = c("phylum" = "phylum_from_phylums"))
  out$kingdom <- out$kingdom_from_phylums
  out$genus <- .fb_coalesce_character(out$species_genus, out$genus_from_genera)
  out$species <- as.character(out$species_epithet)
  out$order <- .fb_coalesce_character(out$order_from_family, out$order_from_orders)
  out$taxonomy_source <- "sealifebase_native"
  out$taxon_name <- ifelse(.fb_non_empty(out$genus) & .fb_non_empty(out$species), paste(out$genus, out$species), NA_character_)

  result <- .order_fishbase_taxonomy_columns(out)
  attr(result, "rank_metadata") <- dplyr::bind_rows(
    .rank_metadata_from_table(gen, "sealifebase", "genus", "genus_from_genera", common_name_col = "genus_common_name"),
    .rank_metadata_from_table(fam, "sealifebase", "family", "family", common_name_col = "family_common_name"),
    .rank_metadata_from_table(ord, "sealifebase", "order", "order_from_orders", common_name_col = "order_common_name"),
    .rank_metadata_from_table(cls, "sealifebase", "class", "class_from_classes", common_name_col = "class_common_name"),
    .rank_metadata_from_table(phy, "sealifebase", "phylum", "phylum_from_phylums", common_name_col = "phylum_common_name")
  )
  result
}

.fb_taxonomy_table_or_empty <- function(taxonomy, source) {
  if (is.null(taxonomy)) {
    return(.empty_fishbase_taxonomy_table(source))
  }
  taxonomy <- tibble::as_tibble(taxonomy)
  expected <- .fishbase_taxonomy_columns()
  for (column_name in expected) {
    if (!(column_name %in% colnames(taxonomy))) {
      taxonomy[[column_name]] <- NA_character_
    }
  }
  taxonomy <- taxonomy[, expected, drop = FALSE]
  rank_metadata <- attr(taxonomy, "rank_metadata", exact = TRUE)
  if (is.null(rank_metadata)) {
    rank_metadata <- .empty_fishbase_rank_metadata()
  }
  attr(taxonomy, "rank_metadata") <- rank_metadata
  taxonomy
}

.empty_fishbase_taxonomy_table <- function(source) {
  out <- tibble::as_tibble(data.frame(
    source = character(),
    spec_code = character(),
    taxon_name_taxonomy = character(),
    kingdom = character(),
    phylum = character(),
    superclass = character(),
    class = character(),
    subclass = character(),
    order = character(),
    family = character(),
    subfamily = character(),
    genus_taxonomy = character(),
    species_taxonomy = character(),
    taxonomy_source = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  attr(out, "rank_metadata") <- .empty_fishbase_rank_metadata()
  out
}

.fishbase_taxonomy_columns <- function() {
  c(
    "source",
    "spec_code",
    "taxon_name_taxonomy",
    "kingdom",
    "phylum",
    "superclass",
    "class",
    "subclass",
    "order",
    "family",
    "subfamily",
    "genus_taxonomy",
    "species_taxonomy",
    "taxonomy_source"
  )
}

.order_fishbase_taxonomy_columns <- function(x) {
  out <- data.frame(
    source = as.character(x$source),
    spec_code = .fb_as_spec_code(x$spec_code),
    taxon_name_taxonomy = as.character(x$taxon_name),
    kingdom = as.character(x$kingdom),
    phylum = as.character(x$phylum),
    superclass = .fb_get_or_na(x, "superclass"),
    class = as.character(x$class),
    subclass = .fb_get_or_na(x, "subclass"),
    order = as.character(x$order),
    family = as.character(x$family),
    subfamily = as.character(x$subfamily),
    genus_taxonomy = as.character(x$genus),
    species_taxonomy = as.character(x$species),
    taxonomy_source = as.character(x$taxonomy_source),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (column_name in colnames(out)) {
    if (is.character(out[[column_name]])) {
      out[[column_name]][!.fb_non_empty(out[[column_name]])] <- NA_character_
    }
  }
  out <- out[!duplicated(paste(out$source, out$spec_code, sep = "\r")), , drop = FALSE]
  tibble::as_tibble(out)
}

.coalesce_joined_taxonomy_columns <- function(x) {
  pairs <- c(
    taxon_name = "taxon_name_taxonomy",
    genus = "genus_taxonomy",
    species = "species_taxonomy"
  )
  for (target in names(pairs)) {
    source_column <- pairs[[target]]
    if (source_column %in% colnames(x)) {
      x[[target]] <- .fb_coalesce_character(x[[source_column]], x[[target]])
      x[[source_column]] <- NULL
    }
  }
  x
}

.rank_metadata_from_table <- function(x, source, taxon_rank, taxon_name_col, common_name_col = NULL) {
  if (is.null(x) || nrow(x) == 0 || !(taxon_name_col %in% colnames(x))) {
    return(.empty_fishbase_rank_metadata())
  }
  taxon_name <- as.character(x[[taxon_name_col]])
  out <- data.frame(
    source = source,
    taxon_rank = taxon_rank,
    taxon_name = taxon_name,
    common_name = .fb_metadata_column(x, common_name_col),
    rank_water_salinity = .fb_metadata_column(x, "water_salinity"),
    rank_body_shape = .fb_metadata_column(x, "body_shape_i"),
    rank_habitat_note = .fb_metadata_column(x, "habitat"),
    rank_distribution = .fb_metadata_column(x, "distribution"),
    rank_diagnosis = .fb_metadata_column(x, "diagnosis"),
    rank_remarks = .collapse_parallel_metadata_columns(x, c("remark", "remarks", "comment", "classification_remark")),
    rank_etymology = .collapse_parallel_metadata_columns(x, .fishbase_etymology_columns()),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[.fb_non_empty(out$taxon_name), , drop = FALSE]
  out <- unique(out)
  tibble::as_tibble(out)
}

.empty_fishbase_rank_metadata <- function() {
  tibble::as_tibble(data.frame(
    source = character(),
    taxon_rank = character(),
    taxon_name = character(),
    common_name = character(),
    rank_water_salinity = character(),
    rank_body_shape = character(),
    rank_habitat_note = character(),
    rank_distribution = character(),
    rank_diagnosis = character(),
    rank_remarks = character(),
    rank_etymology = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}

.fb_metadata_column <- function(x, column_name) {
  if (is.null(column_name) || !(column_name %in% colnames(x))) {
    return(rep(NA_character_, nrow(x)))
  }
  out <- as.character(x[[column_name]])
  out[!.fb_non_empty(out)] <- NA_character_
  out
}

.collapse_parallel_metadata_columns <- function(x, columns) {
  values <- lapply(columns, function(column_name) .fb_metadata_column(x, column_name))
  values <- do.call(cbind, values)
  apply(values, 1, .collapse_fishbase_text_values)
}

.fishbase_etymology_columns <- function() {
  c(
    "etymology",
    "etymology_1",
    "etymology1",
    "etymology_2",
    "etymology2",
    "name_etymology",
    "etym"
  )
}

.order_clean_fishbase_ecology_columns <- function(x) {
  expected <- c(
    "source",
    "spec_code",
    "ecology_comments",
    "ecology_zone",
    "ecology_water_column_zone",
    "ecology_mobility",
    "ecology_size_class",
    "ecology_substrate",
    "ecology_special_habitat"
  )
  for (column_name in expected) {
    if (!(column_name %in% colnames(x))) {
      x[[column_name]] <- NA_character_
    }
  }
  tibble::as_tibble(x[, expected, drop = FALSE])
}

.order_fishbase_ecology_db_columns <- function(x) {
  expected <- .fishbase_species_db_columns()
  for (column_name in expected) {
    if (!(column_name %in% colnames(x))) {
      x[[column_name]] <- NA_character_
    }
  }
  trailing <- setdiff(colnames(x), expected)
  tibble::as_tibble(x[, c(expected, trailing), drop = FALSE])
}

.fishbase_species_db_columns <- function() {
  c(
    "source",
    "spec_code",
    "taxon_rank",
    "taxon_name",
    "kingdom",
    "phylum",
    "superclass",
    "class",
    "subclass",
    "order",
    "family",
    "subfamily",
    "genus",
    "species",
    "common_name",
    "taxonomy_source",
    "ecology_environment",
    "ecology_position",
    "ecology_habitat_broad",
    "ecology_depth_zone",
    "depth_min",
    "depth_max",
    "common_depth_min",
    "common_depth_max",
    "ecology_zone",
    "ecology_water_column_zone",
    "ecology_mobility",
    "ecology_size_class",
    "ecology_substrate",
    "ecology_special_habitat",
    "raw_demers_pelag",
    "species_comments",
    "ecology_comments"
  )
}

.fishbase_environment_text <- function(species) {
  flags <- c(
    saltwater = "marine",
    brack = "brackish",
    fresh = "freshwater",
    land = "terrestrial"
  )
  .collapse_fishbase_flag_category(
    species,
    flags,
    present_fun = .fb_environment_present
  )
}

.fishbase_position_text <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value[value %in% c("", "na", "nan", "null") | is.na(value)] <- NA_character_
  mapped <- dplyr::case_when(
    value %in% c("benthic", "sessile", "reef-associated") ~ "benthic",
    value %in% c("demersal", "bathydemersal") ~ "demersal",
    value == "benthopelagic" ~ "benthopelagic",
    value %in% c("pelagic", "pelagic-neritic", "pelagic-oceanic", "epipelagic", "bathypelagic") ~ "pelagic",
    value %in% c("host", "epiphytic") ~ "host_associated",
    value == "others" ~ "other",
    is.na(value) | value == "unknown" ~ "unknown",
    TRUE ~ "unknown"
  )
  as.character(mapped)
}

.fishbase_depth_zone <- function(depth_max) {
  depth <- suppressWarnings(as.numeric(depth_max))
  out <- rep("unknown", length(depth))
  out[!is.na(depth) & depth < 200] <- "shallow_only"
  out[!is.na(depth) & depth >= 200 & depth < 1000] <- "mesophotic_or_upper_bathyal_possible"
  out[!is.na(depth) & depth >= 1000 & depth < 3000] <- "bathyal_possible"
  out[!is.na(depth) & depth >= 3000 & depth < 6000] <- "abyssal_possible"
  out[!is.na(depth) & depth >= 6000] <- "hadal_or_abyssal_possible"
  out
}

.fishbase_habitat_broad <- function(ecology_position, ecology_habitat_broad) {
  vapply(seq_along(ecology_position), function(index) {
    values <- character()
    position <- as.character(ecology_position[[index]])
    if (position %in% c("benthic", "demersal", "benthopelagic", "pelagic", "host_associated", "other")) {
      values <- c(values, position)
    }
    values <- c(values, .split_fishbase_collapsed_text(ecology_habitat_broad[[index]]))
    values <- unique(values[.fb_non_empty(values) & values != "unknown"])
    allowed <- c("benthic", "demersal", "benthopelagic", "pelagic", "host_associated", "other")
    values <- allowed[allowed %in% values]
    if (length(values) == 0) {
      return("unknown")
    }
    paste(values, collapse = "; ")
  }, character(1))
}

.fishbase_ecology_zone_map <- function() {
  c(
    neritic = "neritic",
    oceanic = "oceanic",
    intertidal = "intertidal",
    sub_littoral = "sublittoral"
  )
}

.fishbase_water_column_zone_map <- function() {
  c(
    epipelagic = "epipelagic",
    mesopelagic = "mesopelagic",
    bathypelagic = "bathypelagic",
    abyssopelagic = "abyssopelagic",
    hadopelagic = "hadopelagic"
  )
}

.fishbase_habitat_broad_map <- function() {
  c(
    benthic = "benthic",
    demersal = "demersal",
    pelagic = "pelagic"
  )
}

.fishbase_mobility_map <- function() {
  c(
    sessile = "sessile",
    mobile = "mobile"
  )
}

.fishbase_size_class_map <- function() {
  c(
    endofauna = "endofauna",
    macrobenthos = "macrobenthos",
    megabenthos = "megabenthos",
    meiobenthos = "meiobenthos"
  )
}

.fishbase_substrate_map <- function() {
  c(
    soft_bottom = "soft_bottom",
    sand = "sand",
    mud = "mud",
    ooze = "ooze",
    hard_bottom = "hard_bottom",
    rocky = "rocky",
    rubble = "rubble",
    gravel = "gravel",
    coral_reefs = "coral_reefs"
  )
}

.fishbase_special_habitat_map <- function() {
  c(
    seamounts = "seamounts",
    cold_seeps = "cold_seeps",
    hydrothermal_vents = "hydrothermal_vents",
    vents = "vents",
    deep_water_corals = "deep_water_corals"
  )
}

.collapse_fishbase_flag_category <- function(x, mapping, present_fun = .fb_flag_present) {
  x <- .fb_add_missing_columns(x, names(mapping))
  present <- vapply(names(mapping), function(column_name) {
    present_fun(x[[column_name]])
  }, logical(nrow(x)))
  if (is.null(dim(present))) {
    present <- matrix(
      present,
      nrow = nrow(x),
      ncol = length(mapping),
      dimnames = list(NULL, names(mapping))
    )
  }
  vapply(seq_len(nrow(x)), function(row_index) {
    values <- unname(mapping[which(present[row_index, ])])
    if (length(values) == 0) {
      return(NA_character_)
    }
    paste(unique(values), collapse = "; ")
  }, character(1))
}

.collapse_fishbase_flag_table <- function(x, source, mapping, output_column) {
  rows <- lapply(names(mapping), function(column_name) {
    present <- .fb_flag_present(x[[column_name]])
    spec_code <- .fb_as_spec_code(x$spec_code)
    use <- present & .fb_non_empty(spec_code)
    if (!any(use)) {
      return(NULL)
    }
    data.frame(
      source = source,
      spec_code = spec_code[use],
      value = unname(mapping[[column_name]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    out <- data.frame(
      source = character(),
      spec_code = character(),
      value = character(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  }
  if (nrow(out) == 0) {
    result <- data.frame(
      source = character(),
      spec_code = character(),
      value = character(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(result)[[3]] <- output_column
    return(result)
  }
  out <- unique(out)
  out <- out[, .(value = paste(sort(unique(value)), collapse = "; ")),
    by = .(source, spec_code)
  ]
  data.table::setnames(out, "value", output_column)
  as.data.frame(out, stringsAsFactors = FALSE)
}

.collapse_fishbase_text_values <- function(values) {
  values <- trimws(as.character(values))
  values <- values[.fb_non_empty(values)]
  if (length(values) == 0) {
    return(NA_character_)
  }
  if (any(grepl(";", values, fixed = TRUE))) {
    split <- unlist(strsplit(values, "\\s*;\\s*"), use.names = FALSE)
  } else {
    split <- values
  }
  split <- trimws(split)
  split <- split[.fb_non_empty(split)]
  split <- sort(unique(split))
  if (length(split) == 0) {
    return(NA_character_)
  }
  paste(split, collapse = "; ")
}

.collapse_fishbase_atomic_values <- function(values) {
  values <- trimws(as.character(values))
  values <- values[.fb_non_empty(values)]
  values <- values[tolower(values) != "unknown"]
  values <- sort(unique(values))
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(values, collapse = "; ")
}

.split_fishbase_collapsed_text <- function(value) {
  value <- as.character(value)
  if (!.fb_non_empty(value)) {
    return(character())
  }
  out <- unlist(strsplit(value, "\\s*;\\s*"), use.names = FALSE)
  trimws(out[.fb_non_empty(out)])
}

.clean_fishbase_names <- function(x, name) {
  if (!inherits(x, "data.frame")) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  colnames(x) <- janitor::make_clean_names(colnames(x))
  x
}

.fb_add_missing_columns <- function(x, columns) {
  for (column_name in columns) {
    if (!(column_name %in% colnames(x))) {
      x[[column_name]] <- NA
    }
  }
  x
}

.fb_null_data_frame <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }
  x
}

.fb_distinct_by_key <- function(x, key) {
  if (!(key %in% colnames(x)) || nrow(x) == 0) {
    return(x)
  }
  x <- x[.fb_non_empty(x[[key]]), , drop = FALSE]
  x[!duplicated(x[[key]]), , drop = FALSE]
}

.fb_get_or_na <- function(x, column_name) {
  if (!(column_name %in% colnames(x))) {
    return(rep(NA_character_, nrow(x)))
  }
  as.character(x[[column_name]])
}

.fb_coalesce_character <- function(...) {
  values <- list(...)
  if (length(values) == 0) {
    return(character())
  }
  max_length <- max(vapply(values, length, integer(1)))
  out <- rep(NA_character_, max_length)
  for (value in values) {
    value <- as.character(value)
    if (length(value) == 1 && max_length > 1) {
      value <- rep(value, max_length)
    }
    value[!.fb_non_empty(value)] <- NA_character_
    use <- !.fb_non_empty(out) & .fb_non_empty(value)
    out[use] <- value[use]
  }
  out
}

.fb_first_non_empty_column <- function(x, columns) {
  out <- rep(NA_character_, nrow(x))
  for (column_name in columns) {
    if (!(column_name %in% colnames(x))) {
      next
    }
    value <- as.character(x[[column_name]])
    use <- !.fb_non_empty(out) & .fb_non_empty(value)
    out[use] <- value[use]
  }
  out
}

.fb_as_spec_code <- function(value) {
  value <- as.character(value)
  value[!.fb_non_empty(value)] <- NA_character_
  value
}

.fb_flag_present <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value %in% c("-1", "1", "true", "t", "yes", "y")
}

.fb_environment_present <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value %in% c("1", "-1", "true", "t", "yes", "y")
}

.fb_non_empty <- function(value) {
  value <- as.character(value)
  !is.na(value) & nzchar(trimws(value)) & !(tolower(trimws(value)) %in% c("na", "nan", "null", "none"))
}

.fb_required_scalar_character <- function(value, name) {
  if (!is.character(value) || length(value) != 1 || is.na(value) || !nzchar(value)) {
    stop("`", name, "` must be a non-empty character scalar.", call. = FALSE)
  }
  value
}

.fb_optional_scalar_character <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }
  .fb_required_scalar_character(value, name)
}

.fb_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.require_arrow_for_fishbase_parquet <- function() {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(
      "Building the FishBase/SeaLifeBase ecology database requires the ",
      "suggested package `arrow`. Install it with install.packages(\"arrow\") ",
      "or read parquet files separately and pass data frames to the cleaning ",
      "helpers for testing.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.require_curl_for_fishbase_download <- function() {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop(
      "Downloading FishBase/SeaLifeBase parquet files requires the suggested ",
      "package `curl`. Install it with install.packages(\"curl\") or download ",
      "the parquet files manually with command-line curl.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
