#' Add taxon-level metadata to a phyloseq taxonomy table
#'
#' @description
#' `add_taxa_metadata()` attaches taxon-level metadata, such as the `summary`
#' table returned by [assess_taxa_evidence()], to the `tax_table()` of a
#' phyloseq object. The metadata are stored as additional taxonomy-table
#' columns, so they will appear in `phyloseq::rank_names(ps)`.
#'
#' A short prefix is added to each metadata column to make these fields easy to
#' identify and remove. Use [drop_taxa_metadata()] before taxonomy filtering or
#' merging steps that should only see true taxonomic ranks.
#'
#' Richer evidence tables such as `details`, `worms`, `obis`, and `ecology` may
#' be better kept separately if they contain multiple rows per `feature_id`.
#' The selected table must contain at most one row per `feature_id`.
#'
#' @param ps A phyloseq object with a `tax_table()`.
#' @param taxa_metadata Either the full list returned by
#'   [assess_taxa_evidence()] or a single metadata data frame.
#' @param source Character scalar naming the list element to use when
#'   `taxa_metadata` is a list. Defaults to `"summary"`.
#' @param cols Optional character vector of metadata columns to add. If `NULL`,
#'   all columns except `feature_id` are added.
#' @param prefix Character scalar prepended to added metadata columns. Defaults
#'   to `"tm_"`.
#' @param overwrite Logical scalar. If `FALSE`, error when a target metadata
#'   column already exists in `phyloseq::tax_table(ps)`. If `TRUE`, replace the
#'   target metadata columns.
#' @param verbose Logical scalar. If `TRUE`, report how many columns were added.
#'
#' @return A standard phyloseq object with metadata columns appended to
#'   `phyloseq::tax_table(ps)`.
#' @export
#'
#' @examples
#' otu <- phyloseq::otu_table(
#'   matrix(
#'     c(1, 0, 2, 3),
#'     nrow = 2,
#'     dimnames = list(c("asv1", "asv2"), c("sample1", "sample2"))
#'   ),
#'   taxa_are_rows = TRUE
#' )
#' tax <- phyloseq::tax_table(
#'   matrix(
#'     c("Metazoa", "Chordata", "Metazoa", "Arthropoda"),
#'     nrow = 2,
#'     byrow = TRUE,
#'     dimnames = list(c("asv1", "asv2"), c("kingdom", "phylum"))
#'   )
#' )
#' ps <- phyloseq::phyloseq(otu, tax)
#' metadata <- data.frame(
#'   feature_id = c("asv1", "asv2"),
#'   recommended_action = c("retain", "flag_for_review")
#' )
#'
#' ps <- add_taxa_metadata(ps, metadata, verbose = FALSE)
#' phyloseq::rank_names(ps)
#' ps <- drop_taxa_metadata(ps)
add_taxa_metadata <- function(
  ps,
  taxa_metadata,
  source = "summary",
  cols = NULL,
  prefix = "tm_",
  overwrite = FALSE,
  verbose = TRUE
) {
  .validate_taxa_metadata_phyloseq(ps)
  source <- .validate_taxa_metadata_scalar_character(source, "source")
  prefix <- .validate_taxa_metadata_prefix(prefix)
  overwrite <- .validate_taxa_metadata_scalar_logical(overwrite, "overwrite")
  verbose <- .validate_taxa_metadata_scalar_logical(verbose, "verbose")

  metadata <- .select_taxa_metadata_table(taxa_metadata, source = source)
  metadata <- .validate_taxa_metadata_table(metadata)
  cols <- .select_taxa_metadata_columns(metadata, cols)

  tax_table <- methods::as(phyloseq::tax_table(ps, errorIfNULL = FALSE), "matrix")
  target_names <- paste0(prefix, cols)
  if (anyDuplicated(target_names)) {
    stop("`cols` produces duplicated metadata column names after prefixing.", call. = FALSE)
  }

  existing_target <- intersect(target_names, colnames(tax_table))
  if (length(existing_target) > 0 && !isTRUE(overwrite)) {
    stop(
      "Metadata column(s) already exist in `tax_table(ps)`: ",
      paste(existing_target, collapse = ", "),
      ". Use `overwrite = TRUE` to replace them.",
      call. = FALSE
    )
  }

  metadata_values <- .build_taxa_metadata_matrix(
    metadata = metadata,
    cols = cols,
    target_names = target_names,
    taxa_ids = phyloseq::taxa_names(ps)
  )

  if (length(existing_target) > 0) {
    tax_table <- tax_table[, !(colnames(tax_table) %in% existing_target), drop = FALSE]
  }

  out_tax_table <- cbind(tax_table, metadata_values)
  rownames(out_tax_table) <- phyloseq::taxa_names(ps)
  phyloseq::tax_table(ps) <- phyloseq::tax_table(out_tax_table)

  if (isTRUE(verbose)) {
    message(
      "Added ",
      length(cols),
      " taxon metadata column",
      if (length(cols) == 1L) "" else "s",
      " to `tax_table(ps)`."
    )
  }

  ps
}

#' Drop taxon-level metadata from a phyloseq taxonomy table
#'
#' @description
#' `drop_taxa_metadata()` removes columns from `phyloseq::tax_table(ps)` whose
#' names start with `prefix`. This is useful because metadata inserted by
#' [add_taxa_metadata()] appears in `phyloseq::rank_names(ps)` and may otherwise
#' be treated as taxonomic ranks by downstream phyloseq workflows.
#'
#' @param ps A phyloseq object with a `tax_table()`.
#' @param prefix Character scalar identifying metadata columns. Defaults to
#'   `"tm_"`.
#'
#' @return A standard phyloseq object with matching metadata columns removed
#'   from `phyloseq::tax_table(ps)`.
#' @export
#'
#' @examples
#' otu <- phyloseq::otu_table(
#'   matrix(
#'     c(1, 0, 2, 3),
#'     nrow = 2,
#'     dimnames = list(c("asv1", "asv2"), c("sample1", "sample2"))
#'   ),
#'   taxa_are_rows = TRUE
#' )
#' tax <- phyloseq::tax_table(
#'   matrix(
#'     c("Metazoa", "Chordata", "Metazoa", "Arthropoda"),
#'     nrow = 2,
#'     byrow = TRUE,
#'     dimnames = list(c("asv1", "asv2"), c("kingdom", "phylum"))
#'   )
#' )
#' ps <- phyloseq::phyloseq(otu, tax)
#' metadata <- data.frame(feature_id = "asv1", recommended_action = "retain")
#' ps <- add_taxa_metadata(ps, metadata, verbose = FALSE)
#' ps <- drop_taxa_metadata(ps)
drop_taxa_metadata <- function(ps, prefix = "tm_") {
  .validate_taxa_metadata_phyloseq(ps)
  prefix <- .validate_taxa_metadata_prefix(prefix)

  tax_table <- methods::as(phyloseq::tax_table(ps, errorIfNULL = FALSE), "matrix")
  keep <- !startsWith(colnames(tax_table), prefix)
  if (all(keep)) {
    return(ps)
  }

  out_tax_table <- tax_table[, keep, drop = FALSE]
  rownames(out_tax_table) <- phyloseq::taxa_names(ps)
  phyloseq::tax_table(ps) <- phyloseq::tax_table(out_tax_table)
  ps
}

.validate_taxa_metadata_phyloseq <- function(ps) {
  if (!methods::is(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }
  if (is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    stop("`ps` must contain a `tax_table()`.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_taxa_metadata_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a single non-empty character string.", call. = FALSE)
  }
  x
}

.validate_taxa_metadata_prefix <- function(prefix) {
  .validate_taxa_metadata_scalar_character(prefix, "prefix")
}

.validate_taxa_metadata_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.select_taxa_metadata_table <- function(taxa_metadata, source) {
  if (inherits(taxa_metadata, "data.frame")) {
    return(as.data.frame(taxa_metadata, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (!is.list(taxa_metadata)) {
    stop("`taxa_metadata` must be a data frame or a list of data frames.", call. = FALSE)
  }
  if (!(source %in% names(taxa_metadata))) {
    stop("`taxa_metadata` does not contain source `", source, "`.", call. = FALSE)
  }
  selected <- taxa_metadata[[source]]
  if (!inherits(selected, "data.frame")) {
    stop("`taxa_metadata[[source]]` must be a data frame.", call. = FALSE)
  }
  as.data.frame(selected, stringsAsFactors = FALSE, check.names = FALSE)
}

.validate_taxa_metadata_table <- function(metadata) {
  if (!("feature_id" %in% colnames(metadata))) {
    stop("The selected metadata table must contain `feature_id`.", call. = FALSE)
  }
  feature_id <- as.character(metadata$feature_id)
  if (any(is.na(feature_id) | !nzchar(feature_id))) {
    stop("`feature_id` values must be non-missing.", call. = FALSE)
  }
  duplicated_ids <- unique(feature_id[duplicated(feature_id)])
  if (length(duplicated_ids) > 0) {
    stop(
      "Duplicate `feature_id` value(s) found: ",
      paste(duplicated_ids, collapse = ", "),
      ". The selected metadata table must have at most one row per `feature_id`.",
      call. = FALSE
    )
  }
  metadata$feature_id <- feature_id
  metadata
}

.select_taxa_metadata_columns <- function(metadata, cols) {
  if (is.null(cols)) {
    cols <- setdiff(colnames(metadata), "feature_id")
  } else if (!is.character(cols) || any(is.na(cols)) || any(!nzchar(cols))) {
    stop("`cols` must be a character vector of column names.", call. = FALSE)
  }
  if ("feature_id" %in% cols) {
    stop("`cols` should not include `feature_id`.", call. = FALSE)
  }
  if (anyDuplicated(cols)) {
    stop("`cols` contains duplicated column names.", call. = FALSE)
  }
  missing_cols <- setdiff(cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(
      "Requested metadata column(s) not found: ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (length(cols) == 0) {
    stop("No metadata columns were selected.", call. = FALSE)
  }
  cols
}

.build_taxa_metadata_matrix <- function(metadata, cols, target_names, taxa_ids) {
  metadata_ids <- as.character(metadata$feature_id)
  row_index <- match(taxa_ids, metadata_ids)
  out <- matrix(
    NA_character_,
    nrow = length(taxa_ids),
    ncol = length(cols),
    dimnames = list(taxa_ids, target_names)
  )

  matched <- !is.na(row_index)
  for (i in seq_along(cols)) {
    values <- .taxa_metadata_column_as_character(metadata[[cols[[i]]]])
    out[matched, i] <- values[row_index[matched]]
  }

  out
}

.taxa_metadata_column_as_character <- function(x) {
  if (is.list(x) && !inherits(x, "data.frame")) {
    return(vapply(x, .collapse_taxa_metadata_value, character(1)))
  }
  as.character(x)
}

.collapse_taxa_metadata_value <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }
  if (is.data.frame(value) || is.list(value)) {
    value <- unlist(value, recursive = TRUE, use.names = FALSE)
  }
  value <- as.character(value)
  value <- value[!is.na(value)]
  if (length(value) == 0) {
    return(NA_character_)
  }
  paste(value, collapse = "; ")
}
