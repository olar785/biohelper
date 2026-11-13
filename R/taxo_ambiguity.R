.biohelper_resolve_ambiguity <- function(matches, candidates, rows, action, sqlFile) {
  diagnostics <- data.frame(
    feature_id = matches$feature_id[rows],
    candidate_ids = vapply(candidates[rows], paste, character(1), collapse = ","),
    action = rep(action, length(rows)), selected_id = matches$id[rows],
    lca_id = rep(NA_character_, length(rows)), lca_name = rep(NA_character_, length(rows)),
    lca_rank = rep(NA_character_, length(rows)), stringsAsFactors = FALSE
  )
  ids <- matches$id
  if (!length(rows)) return(list(ids = ids, diagnostics = diagnostics))
  if (action == "unassigned") ids[rows] <- NA_character_
  if (action == "lca") {
    query_ids <- unique(unlist(candidates[rows]))
    paths <- .biohelper_ncbi_id_paths(query_ids, sqlFile)
    if (!all(vapply(paths, function(path) path$complete, logical(1)))) {
      stop("Cannot calculate an LCA: ancestry is unavailable for one or more candidate TaxIDs.", call. = FALSE)
    }
    lcas <- vapply(candidates[rows], function(candidate) {
      # Paths run from the candidate towards the root. Their first shared ID
      # is the deepest common ancestor, including intermediate/unranked nodes.
      common <- Reduce(intersect, lapply(candidate, function(id) paths[[id]]$ids))
      paths[[candidate[[1]]]]$ids[paths[[candidate[[1]]]]$ids %in% common][[1]]
    }, character(1))
    unique_lcas <- unique(lcas)
    labels <- stats::setNames(taxonomizr::getRawTaxonomy(unique_lcas, sqlFile), unique_lcas)
    diagnostics$lca_id <- lcas
    diagnostics$lca_name <- vapply(lcas, function(id) {
      if (id == "1") "root" else unname(labels[[id]][[1]])
    }, character(1))
    diagnostics$lca_rank <- vapply(lcas, function(id) {
      if (id == "1") "no rank" else names(labels[[id]])[[1]]
    }, character(1))
    ids[rows] <- lcas
  }
  diagnostics$selected_id <- ids[rows]
  list(ids = ids, diagnostics = diagnostics)
}
