.biohelper_normalisation_cache <- function() {
  list(names = new.env(parent = emptyenv()), lineages = new.env(parent = emptyenv()))
}

.biohelper_taxonomy_profile_keys <- function(df) {
  fields <- lapply(df, function(x) {
    x <- as.character(x)
    ifelse(is.na(x), "N", paste0(nchar(x, type = "bytes"), ":", x))
  })
  do.call(paste, c(fields, list(sep = "|")))
}

.biohelper_cached_ncbi_names <- function(taxa, sqlFile, cache) {
  missing <- unique(taxa[!vapply(taxa, exists, logical(1), envir = cache$names, inherits = FALSE)])
  if (length(missing)) {
    ids <- as.character(withCallingHandlers(
      taxonomizr::getId(taxa = missing, sqlFile = sqlFile, onlyScientific = TRUE),
      warning = function(w) {
        # Multiple name matches are expected here. Report ambiguity after
        # ordered lineage matching and optional BLAST have run.
        message <- conditionMessage(w)
        if (startsWith(message, "Multiple taxa ids found for ") &&
            endsWith(message, ". Collapsing with commas")) {
          invokeRestart("muffleWarning")
        }
      }
    ))
    for (i in seq_along(missing)) assign(missing[[i]], ids[[i]], envir = cache$names)
  }
  unname(vapply(taxa, get, character(1), envir = cache$names, inherits = FALSE))
}

.biohelper_cached_ncbi_lineages <- function(ids, sqlFile, cache) {
  missing <- ids[!vapply(ids, exists, logical(1), envir = cache$lineages, inherits = FALSE)]
  if (length(missing)) {
    lineages <- taxonomizr::getRawTaxonomy(missing, sqlFile = sqlFile)
    for (i in seq_along(missing)) assign(missing[[i]], lineages[[i]], envir = cache$lineages)
  }
  stats::setNames(lapply(ids, get, envir = cache$lineages, inherits = FALSE), ids)
}

.biohelper_blast_message <- function(control, ...) {
  if (isTRUE(control$verbose)) message(...)
}

.biohelper_blast_target_rows <- function(matches, candidates, mode, feature_ids = NULL) {
  rows <- which(lengths(candidates) > if (mode == "ambiguous") 1L else 0L)
  if (!is.null(feature_ids)) {
    unknown <- setdiff(feature_ids, matches$feature_id)
    if (length(unknown)) stop("Unknown BLAST feature IDs: ", paste(utils::head(unknown, 3), collapse = ", "), call. = FALSE)
    rows <- rows[matches$feature_id[rows] %in% feature_ids]
  }
  rows
}

.biohelper_pool_normalisation_blast <- function(plans, sequences, control, db, blastapp_path, sqlFile) {
  all_ids <- unique(unlist(lapply(plans, function(x) x$res_df$feature_id)))
  unknown <- setdiff(control$feature_ids, all_ids)
  if (length(unknown)) stop("Unknown BLAST feature IDs: ", paste(utils::head(unknown, 3), collapse = ", "), call. = FALSE)
  target_rows <- lapply(plans, function(x) .biohelper_blast_target_rows(
    x$res_df, x$blast_candidates, control$scope,
    if (is.null(control$feature_ids)) NULL else intersect(control$feature_ids, x$res_df$feature_id)
  ))
  targets <- unique(unlist(lapply(seq_along(plans), function(i) {
    plans[[i]]$res_df$feature_id[target_rows[[i]]]
  })))
  ambiguous_ids <- unique(unlist(lapply(plans, function(x) {
    x$res_df$feature_id[lengths(x$blast_candidates) > 1L]
  })))
  .biohelper_blast_message(control, "Ambiguous ASVs found: ", length(ambiguous_ids),
    " distinct across all input tables.")
  .biohelper_blast_message(control, "BLAST planning: ", sum(lengths(target_rows)),
    " selected assignments across ", length(plans), " tables; ", length(targets), " distinct ASVs.")
  available <- if (!length(targets) || is.null(sequences)) stats::setNames(character(), character()) else {
    .biohelper_blast_message(control, "Reading ASV sequences once for all input tables.")
    .biohelper_read_blast_sequences(sequences)
  }
  missing <- setdiff(targets, names(available))
  searchable <- intersect(targets, names(available))
  .biohelper_blast_message(control, length(missing), " selected ASVs have no sequence; fallback: '", control$fallback, "'.")
  control$db <- if (is.null(db) && control$remote) "core_nt" else db
  control$blastapp_path <- blastapp_path
  hits <- if (length(searchable)) .biohelper_blast_search(
    .biohelper_blast_sequences(available, searchable), control
  ) else .biohelper_empty_blast_hits()
  qualified <- .biohelper_qualified_blast_hits(hits, control)
  ids <- unique(c(unlist(lapply(seq_along(plans), function(i) {
    rows <- target_rows[[i]]
    if (control$fallback == "lca") rows <- union(rows, which(lengths(plans[[i]]$blast_candidates) > 1L))
    unlist(plans[[i]]$blast_candidates[rows])
  })), trimws(unlist(strsplit(qualified$staxids, "[;,]")))))
  ids <- ids[!is.na(ids) & grepl("^[1-9][0-9]*$", ids)]
  .biohelper_blast_message(control, "Loading NCBI ancestry once for ", length(ids), " candidate/hit TaxIDs.")
  paths <- if (length(ids)) .biohelper_ncbi_id_paths(ids, sqlFile) else list()
  list(hits = hits, missing = missing, paths = paths)
}
