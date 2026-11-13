.biohelper_blast_options <- function(control) {
  defaults <- list(
    remote = FALSE,
    task = "blastn", cache_dir = file.path(tempdir(), "biohelper-blast"), refresh = FALSE,
    max_hits = 100L, evalue = 1e-10, min_identity = 95, min_coverage = 90,
    score_fraction = 0.98, min_accessions = 2L, min_taxa = 2L, min_support = 1,
    max_remote_queries = 25L, nthreads = 1L, feature_ids = NULL, conflict_action = "flag",
    poll_seconds = 60, timeout_seconds = 900, email = NULL,
    scope = "ambiguous", fallback = "first"
  )
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control))))) {
    stop("`blast_options` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown BLAST control(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  for (name in names(control)) defaults[name] <- control[name]
  for (name in c("remote", "refresh")) {
    x <- defaults[[name]]
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("Invalid BLAST control: ", name, call. = FALSE)
  }
  for (name in c("task", "cache_dir", "conflict_action", "scope", "fallback")) {
    x <- defaults[[name]]
    if (!is.character(x) || length(x) != 1L || is.na(x)) stop("Invalid BLAST control: ", name, call. = FALSE)
  }
  if (!nzchar(defaults$cache_dir) ||
      !defaults$task %in% c("blastn", "megablast") ||
      !defaults$conflict_action %in% c("flag", "unassign") ||
      !defaults$scope %in% c("ambiguous", "all") ||
      !defaults$fallback %in% c("first", "lca", "unassigned")) {
    stop("Invalid BLAST task, cache directory, conflict action, scope or fallback.", call. = FALSE)
  }
  for (name in c("max_hits", "min_accessions", "min_taxa", "max_remote_queries", "nthreads")) {
    x <- defaults[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)) {
      stop("BLAST control `", name, "` must be a positive integer.", call. = FALSE)
    }
  }
  for (name in c("min_identity", "min_coverage", "score_fraction", "min_support", "evalue")) {
    x <- defaults[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop("Invalid BLAST control: ", name, call. = FALSE)
    }
  }
  if (defaults$min_identity > 100 || defaults$min_coverage > 100 ||
      defaults$score_fraction > 1 || defaults$min_support > 1 || defaults$min_support <= 0.5) {
    stop("Identity/coverage must be <=100; score_fraction <=1; min_support must be >0.5 and <=1.", call. = FALSE)
  }
  if (!is.null(defaults$feature_ids) && (!is.character(defaults$feature_ids) ||
      anyNA(defaults$feature_ids) || any(!nzchar(defaults$feature_ids)))) {
    stop("`blast_options$feature_ids` must contain feature identifiers.", call. = FALSE)
  }
  if (!is.numeric(defaults$poll_seconds) || length(defaults$poll_seconds) != 1L ||
      !is.finite(defaults$poll_seconds) || defaults$poll_seconds < 60 ||
      !is.numeric(defaults$timeout_seconds) || length(defaults$timeout_seconds) != 1L ||
      !is.finite(defaults$timeout_seconds) || defaults$timeout_seconds <= 0 ||
      (!is.null(defaults$email) && (!is.character(defaults$email) || length(defaults$email) != 1L || is.na(defaults$email)))) {
    stop("BLAST polling must be at least 60 seconds; timeout must be positive; email must be a string.", call. = FALSE)
  }
  defaults
}

.biohelper_read_blast_sequences <- function(sequences) {
  if (is.character(sequences) && length(sequences) == 1L && is.null(names(sequences)) && file.exists(sequences)) {
    sequences <- Biostrings::readDNAStringSet(sequences)
    names(sequences) <- sub("\\s.*$", "", names(sequences))
  }
  if (inherits(sequences, "DNAStringSet")) {
    sequence_names <- names(sequences)
    sequences <- stats::setNames(as.character(sequences), sequence_names)
  }
  if (!is.character(sequences) || is.null(names(sequences)) || anyNA(names(sequences)) ||
      any(!nzchar(names(sequences))) || anyDuplicated(names(sequences))) {
    stop("Supply `sequences` as a FASTA path, named DNAStringSet or uniquely named character vector.", call. = FALSE)
  }
  sequences
}

.biohelper_blast_sequences <- function(sequences, feature_ids) {
  sequences <- .biohelper_read_blast_sequences(sequences)
  missing <- setdiff(feature_ids, names(sequences))
  if (length(missing)) {
    stop("Sequences are missing for ", length(missing), " selected feature(s), e.g. ",
         paste(utils::head(missing, 3L), collapse = ", "),
         ". Restrict blast_options$feature_ids or supply the matching sequences.", call. = FALSE)
  }
  sequences <- toupper(sequences[feature_ids])
  if (anyNA(sequences) || any(!grepl("^[ACGTRYSWKMBDHVN]+$", sequences))) {
    stop("Selected sequences must be non-empty, ungapped IUPAC DNA sequences.", call. = FALSE)
  }
  sequences
}

.biohelper_empty_blast_hits <- function() {
  data.frame(qseqid = character(), sseqid = character(), staxids = character(),
             pident = numeric(), qcovhsp = numeric(), evalue = numeric(), bitscore = numeric())
}

.biohelper_validate_blast_hits <- function(hits) {
  if (is.character(hits) && length(hits) == 1L) hits <- read_blast_taxo_csv(hits)
  required <- names(.biohelper_empty_blast_hits())
  if (!is.data.frame(hits) || !all(required %in% names(hits))) {
    stop("`blast_hits` must be a raw BLAST table with columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  hits <- as.data.frame(hits, stringsAsFactors = FALSE)
  for (name in c("qseqid", "sseqid", "staxids")) hits[[name]] <- as.character(hits[[name]])
  if (anyNA(hits$qseqid) || any(!nzchar(hits$qseqid)) || anyNA(hits$sseqid) || any(!nzchar(hits$sseqid))) {
    stop("BLAST query and subject identifiers must be non-missing.", call. = FALSE)
  }
  for (name in c("pident", "qcovhsp", "evalue", "bitscore")) {
    value <- suppressWarnings(as.numeric(hits[[name]]))
    if (any(!is.finite(value)) || any(value < 0) ||
        (name %in% c("pident", "qcovhsp") && any(value > 100))) {
      stop("Invalid numeric values in BLAST column `", name, "`.", call. = FALSE)
    }
    hits[[name]] <- value
  }
  hits
}

.biohelper_blast_key <- function(sequence, specification) {
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(list(sequence = sequence, specification = specification), path, version = 2)
  unname(tools::md5sum(path))
}

.biohelper_run_blast <- function(query, output, log, control) {
  if (!is.character(control$db) || length(control$db) != 1L || is.na(control$db) || !nzchar(control$db)) {
    stop("Supply `db` as the local BLAST database prefix (for example '/databases/nt'). This is separate from the taxonomy SQLite file in `sqlFile`.", call. = FALSE)
  }
  if (control$remote) return(.biohelper_run_remote_blast(query, output, log, control))
  if (!is.character(control$blastapp_path) || length(control$blastapp_path) != 1L ||
      is.na(control$blastapp_path) || !nzchar(control$blastapp_path) || !file.exists(control$blastapp_path)) {
    stop("NCBI BLAST+ `blastn` was not found. Supply `blastapp_path` or precomputed `blast_hits`.", call. = FALSE)
  }
  output_format <- "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"
  args <- c("-query", shQuote(query), "-db", shQuote(control$db),
            "-task", control$task, "-evalue", format(control$evalue, scientific = TRUE),
            "-max_target_seqs", control$max_hits, "-outfmt", shQuote(output_format), "-out", shQuote(output))
  args <- c(args, "-num_threads", control$nthreads)
  status <- system2(control$blastapp_path, args, stdout = log, stderr = log)
  if (status != 0L || !file.exists(output)) {
    stop("BLAST failed. See ", log, ". No taxonomy has been changed by this search.", call. = FALSE)
  }
  invisible(output)
}

.biohelper_blast_search <- function(sequences, control) {
  # Preserve cache keys from the previous working API, where db was database.
  specification <- c(list(database = control$db), control[c("remote", "task", "max_hits", "evalue")])
  keys <- vapply(unname(sequences), .biohelper_blast_key, character(1), specification = specification)
  unique_keys <- unique(keys)
  unique_sequences <- stats::setNames(unname(sequences[match(unique_keys, keys)]), unique_keys)
  dir.create(control$cache_dir, recursive = TRUE, showWarnings = FALSE)
  cached <- stats::setNames(vector("list", length(unique_keys)), unique_keys)
  for (key in unique_keys) {
    file <- file.path(control$cache_dir, paste0(key, ".rds"))
    if (!control$refresh && file.exists(file)) {
      entry <- readRDS(file)
      if (!identical(entry$sequence, unique_sequences[[key]]) || !identical(entry$specification, specification)) {
        stop("BLAST cache metadata does not match its key: ", file, call. = FALSE)
      }
      cached[[key]] <- entry$hits
    }
  }
  pending <- names(cached)[vapply(cached, is.null, logical(1))]
  if (control$remote && length(pending) > control$max_remote_queries) {
    stop("This needs ", length(pending), " uncached remote queries; the configured limit is ",
         control$max_remote_queries, ". Select a pilot with blast_options$feature_ids, use local BLAST, or supply precomputed blast_hits.", call. = FALSE)
  }
  if (length(pending)) {
    batch_key <- .biohelper_blast_key(paste(sort(pending), collapse = ":"), specification)
    prefix <- file.path(control$cache_dir, paste0("batch_", batch_key))
    query <- paste0(prefix, ".fasta")
    output <- paste0(prefix, ".tsv")
    dna <- Biostrings::DNAStringSet(unique_sequences[pending])
    names(dna) <- paste0("bh_", pending)
    Biostrings::writeXStringSet(dna, query)
    .biohelper_run_blast(query, output, paste0(prefix, ".log"), control)
    hits <- if (file.info(output)$size == 0) .biohelper_empty_blast_hits() else .biohelper_validate_blast_hits(output)
    for (key in pending) {
      cached[[key]] <- hits[hits$qseqid == paste0("bh_", key), , drop = FALSE]
      saveRDS(list(sequence = unique_sequences[[key]], specification = specification,
                   retrieved_utc = format(Sys.time(), tz = "UTC"), hits = cached[[key]]),
              file.path(control$cache_dir, paste0(key, ".rds")))
    }
  }
  output <- lapply(seq_along(sequences), function(i) {
    hits <- cached[[keys[[i]]]]
    hits$qseqid <- rep(names(sequences)[[i]], nrow(hits))
    hits
  })
  if (length(output)) do.call(rbind, output) else .biohelper_empty_blast_hits()
}

# TaxIDs, rather than names or rank labels, identify membership in a lineage.
.biohelper_ncbi_id_paths <- function(ids, sqlFile) {
  ids <- unique(as.character(ids))
  connection <- DBI::dbConnect(RSQLite::SQLite(), sqlFile, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(connection))
  parents <- new.env(parent = emptyenv())
  assign("1", "1", parents)
  pending <- setdiff(ids, "1")
  depth <- 0L
  while (length(pending)) {
    next_ids <- character()
    batches <- split(pending, ceiling(seq_along(pending) / 900L))
    for (batch in batches) {
      sql <- paste0("SELECT id, parent FROM nodes WHERE id IN (", paste(rep("?", length(batch)), collapse = ","), ")")
      rows <- DBI::dbGetQuery(connection, sql, params = as.list(as.numeric(batch)))
      found <- as.character(rows$id)
      for (j in seq_along(batch)) {
        index <- match(batch[[j]], found)
        parent <- if (is.na(index)) NA_character_ else as.character(rows$parent[[index]])
        assign(batch[[j]], parent, parents)
        if (!is.na(parent)) next_ids <- c(next_ids, parent)
      }
    }
    pending <- unique(next_ids[!vapply(next_ids, exists, logical(1), envir = parents, inherits = FALSE)])
    depth <- depth + 1L
    if (depth > 200L) stop("NCBI taxonomy ancestry exceeded 200 levels.", call. = FALSE)
  }
  stats::setNames(lapply(ids, function(id) {
    path <- character()
    while (!is.na(id) && !id %in% path) {
      path <- c(path, id)
      if (id == "1") return(list(ids = path, complete = TRUE))
      id <- get(id, parents, inherits = FALSE)
    }
    list(ids = path, complete = FALSE)
  }), ids)
}

.biohelper_check_blast <- function(matches, candidates, sqlFile, mode, sequences, hits, control,
                                   db = NULL, blastapp_path = unname(Sys.which("blastn"))) {
  hit_source <- if (is.null(hits)) "blastn" else "precomputed"
  control <- .biohelper_blast_options(control)
  control$db <- if (is.null(db) && control$remote) "core_nt" else db
  control$blastapp_path <- blastapp_path
  rows <- which(lengths(candidates) > if (mode == "ambiguous") 1L else 0L)
  if (!is.null(control$feature_ids)) {
    unknown <- setdiff(control$feature_ids, matches$feature_id)
    if (length(unknown)) stop("Unknown BLAST feature IDs: ", paste(utils::head(unknown, 3), collapse = ", "), call. = FALSE)
    rows <- rows[matches$feature_id[rows] %in% control$feature_ids]
  }
  target <- matches$feature_id[rows]
  missing_sequences <- character()
  if (is.null(hits)) {
    if (length(rows)) {
      available <- if (is.null(sequences)) stats::setNames(character(), character()) else .biohelper_read_blast_sequences(sequences)
      missing_sequences <- setdiff(target, names(available))
      searchable <- intersect(target, names(available))
      hits <- if (length(searchable)) .biohelper_blast_search(
        .biohelper_blast_sequences(available, searchable), control
      ) else .biohelper_empty_blast_hits()
    } else hits <- .biohelper_empty_blast_hits()
  } else hits <- .biohelper_validate_blast_hits(hits)
  hits <- hits[hits$qseqid %in% target, , drop = FALSE]
  quality <- hits$pident >= control$min_identity & hits$qcovhsp >= control$min_coverage & hits$evalue <= control$evalue
  eligible <- vector("list", length(rows))
  taxids <- vector("list", length(rows))
  for (j in seq_along(rows)) {
    x <- hits[hits$qseqid == target[[j]] & quality, , drop = FALSE]
    if (nrow(x)) {
      x <- x[x$bitscore >= max(x$bitscore) * control$score_fraction, , drop = FALSE]
      x <- x[order(-x$bitscore), , drop = FALSE]
      x <- x[!duplicated(x$sseqid), , drop = FALSE]
    }
    eligible[[j]] <- x
    taxids[[j]] <- unique(trimws(unlist(strsplit(x$staxids, "[;,]"))))
  }
  valid_ids <- unique(c(unlist(candidates[rows]), unlist(taxids)))
  valid_ids <- valid_ids[!is.na(valid_ids) & grepl("^[1-9][0-9]*$", valid_ids)]
  paths <- .biohelper_ncbi_id_paths(valid_ids, sqlFile)
  ids <- matches$id
  diagnostics <- vector("list", length(rows))
  for (j in seq_along(rows)) {
    row <- rows[[j]]
    candidate <- candidates[[row]]
    x <- eligible[[j]]
    taxa <- taxids[[j]]
    selected <- NA_character_
    support <- NA_real_
    status <- if (target[[j]] %in% missing_sequences) "missing_sequence" else "no_hits"
    if (nrow(x)) {
      available <- vapply(c(candidate, taxa), function(id) {
        !is.na(id) && id %in% names(paths) && paths[[id]]$complete
      }, logical(1))
      if (!all(available) || !length(taxa)) status <- "taxonomy_unavailable"
      else if (nrow(x) < control$min_accessions || length(taxa) < control$min_taxa) status <- "insufficient_records"
      else {
        fractions <- vapply(candidate, function(id) mean(vapply(taxa, function(taxon) id %in% paths[[taxon]]$ids, logical(1))), numeric(1))
        support <- max(fractions)
        winners <- candidate[fractions >= control$min_support]
        if (length(winners)) {
          deepest <- winners[which.max(vapply(winners, function(id) length(paths[[id]]$ids), integer(1)))]
          if (all(winners %in% paths[[deepest]]$ids)) {
            selected <- deepest
            status <- if (length(candidate) > 1L) "blast_resolved" else "blast_supported"
            ids[[row]] <- selected
          } else status <- "mixed_hits"
        } else {
          compatible <- vapply(taxa, function(taxon) any(vapply(candidate, function(id) {
            id %in% paths[[taxon]]$ids || taxon %in% paths[[id]]$ids
          }, logical(1))), logical(1))
          status <- if (!any(compatible)) "blast_conflict" else if (all(compatible) && support == 0) "insufficient_resolution" else "mixed_hits"
          if (status == "blast_conflict" && control$conflict_action == "unassign") ids[[row]] <- NA_character_
        }
      }
    } else if (any(hits$qseqid == target[[j]])) status <- "insufficient_quality"
    diagnostics[[j]] <- data.frame(
      feature_id = target[[j]], candidate_ids = paste(candidate, collapse = ","),
      name_selected_id = matches$id[[row]], blast_selected_id = selected, final_id = ids[[row]],
      status = status, n_raw_hits = sum(hits$qseqid == target[[j]]), n_accessions = nrow(x),
      n_taxa = length(taxa), support_fraction = support,
      min_identity = if (nrow(x)) min(x$pident) else NA_real_,
      min_coverage = if (nrow(x)) min(x$qcovhsp) else NA_real_, stringsAsFactors = FALSE
    )
  }
  diagnostics <- if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame()
  list(ids = ids, diagnostics = diagnostics, hits = hits, settings = control, hit_source = hit_source)
}
