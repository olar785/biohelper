.biohelper_blast_stats_columns <- function() {
  c("pident", "qcovs", "evalue", "assignment_method", "support_method")
}

.biohelper_blast_raw_columns <- function() {
  c("qseqid", "qlen", "pident", "sseqid", "sgi", "evalue", "bitscore",
    "staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp")
}

.biohelper_read_blast_stage <- function(path, allow_empty = TRUE) {
  if (!file.exists(path)) stop("BLAST output is missing: ", path, call. = FALSE)
  if (file.info(path)$size == 0) {
    if (!allow_empty) stop("An empty legacy BLAST file cannot establish successful completion.", call. = FALSE)
    columns <- .biohelper_blast_raw_columns()
    out <- as.data.frame(stats::setNames(rep(list(character()), length(columns)), columns))
    for (column in c("qlen", "pident", "evalue", "bitscore", "qcovs", "qcovhsp")) out[[column]] <- numeric()
    return(out)
  }
  as.data.frame(read_blast_taxo_csv(path), stringsAsFactors = FALSE)
}

.biohelper_representative_order <- function(x, method = NULL) {
  number <- function(name, decreasing = FALSE) {
    value <- if (name %in% names(x)) suppressWarnings(as.numeric(x[[name]])) else rep(NA_real_, nrow(x))
    value[is.na(value)] <- if (decreasing) -Inf else Inf
    if (decreasing) -value else value
  }
  subject <- if ("sseqid" %in% names(x)) as.character(x$sseqid) else rep("", nrow(x))
  subject_missing <- is.na(subject) | !nzchar(subject)
  subject[is.na(subject)] <- ""
  # Remaining ties are resolved by a canonical row representation, never input order.
  canonical <- vapply(seq_len(nrow(x)), function(i) paste(as.character(x[i, sort(names(x)), drop = TRUE]), collapse = "\034"), character(1))
  method_order <- if (is.null(method)) rep(0L, nrow(x)) else match(method, c("megablast", "blastn"))
  order(number("evalue"), number("bitscore", TRUE), number("pident", TRUE),
        number("qcovs", TRUE), subject_missing, subject, method_order, canonical, na.last = TRUE, method = "radix")
}

.biohelper_file_fingerprint <- function(path, content = TRUE) {
  path <- path.expand(path)
  info <- file.info(path)
  if (!file.exists(path) || isTRUE(info$isdir)) return(list(path = normalizePath(path, mustWork = FALSE), exists = FALSE))
  list(path = normalizePath(path, mustWork = TRUE), exists = TRUE, size = unname(info$size),
       mtime = if (content) NULL else as.numeric(info$mtime), md5 = if (content) unname(tools::md5sum(path)) else NULL)
}

.biohelper_object_hash <- function(x) {
  file <- tempfile()
  on.exit(unlink(file))
  saveRDS(x, file, version = 2)
  unname(tools::md5sum(file))
}

.biohelper_implementation_hash <- function(functions) {
  # Source-reference environments/bytecode are not stable across R sessions.
  # Hash only the actual implementation text, excluding those attributes.
  .biohelper_object_hash(lapply(functions, function(fun) {
    deparse(body(fun), width.cutoff = 500L, control = c("keepNA", "keepInteger"))
  }))
}

.biohelper_acquire_workflow_lock <- function(directory) {
  parent <- file.path(directory, ".biohelper")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(parent, "run.lock")
  owner_file <- file.path(lock, "owner.rds")
  host <- unname(Sys.info()[["nodename"]])
  read_owner <- function() tryCatch(readRDS(owner_file), error = function(e) NULL)
  if (!dir.create(lock, showWarnings = FALSE)) {
    owner <- read_owner()
    same_host <- !is.null(owner) && identical(owner$host, host)
    alive <- if (same_host) tryCatch(isTRUE(suppressWarnings(tools::pskill(owner$pid, signal = 0L))), error = function(e) TRUE) else TRUE
    if (!same_host || alive) stop("Output directory is locked by an active or unverifiable workflow. Use a separate output directory; remove ", lock, " only after confirming the other process has stopped.", call. = FALSE)
    unlink(lock, recursive = TRUE)
    if (!dir.create(lock, showWarnings = FALSE)) stop("Cannot recover the interrupted workflow lock.", call. = FALSE)
    message("Recovered an interrupted workflow lock.")
  }
  token <- .biohelper_object_hash(list(Sys.getpid(), as.numeric(Sys.time()), tempfile()))
  tryCatch(saveRDS(list(pid = Sys.getpid(), host = host, token = token), owner_file),
    error = function(e) { unlink(lock, recursive = TRUE); stop(e) })
  function() {
    owner <- read_owner()
    if (!is.null(owner) && identical(owner$token, token)) unlink(lock, recursive = TRUE)
    invisible(NULL)
  }
}

.biohelper_atomic_publish <- function(from, to, overwrite = FALSE, archive_dir = NULL, provenance = NULL) {
  if (!file.exists(from)) stop("Completed temporary output is missing: ", from, call. = FALSE)
  backup <- NULL
  linked <- FALSE
  if (file.exists(to)) {
    # POSIX rename can replace atomically. For preserved results, a hard link
    # retains the old inode without copying a potentially large BLAST file.
    if (overwrite && isTRUE(suppressWarnings(file.rename(from, to)))) return(invisible(to))
    if (is.null(archive_dir)) archive_dir <- file.path(dirname(to), ".biohelper", "archive")
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    backup <- tempfile(paste0(basename(to), "."), tmpdir = archive_dir)
    linked <- isTRUE(suppressWarnings(file.link(to, backup)))
    if (!linked && !file.rename(to, backup)) stop("Cannot preserve existing output: ", to, call. = FALSE)
  }
  if (!isTRUE(suppressWarnings(file.rename(from, to)))) {
    if (!is.null(backup)) {
      if (linked) unlink(backup) else file.rename(backup, to)
    }
    stop("Cannot publish completed output: ", to, call. = FALSE)
  }
  if (!is.null(backup)) {
    if (overwrite) unlink(backup) else {
      if (!is.null(provenance)) saveRDS(provenance, paste0(backup, ".checkpoint.rds"), version = 2)
      message("Preserved previous file: ", backup)
    }
  }
  invisible(to)
}

.biohelper_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("metadata_", tmpdir = dirname(path))
  on.exit(unlink(temporary))
  saveRDS(x, temporary, version = 2)
  # Metadata is agent-managed; the existing data files are preserved separately.
  .biohelper_atomic_publish(temporary, path, overwrite = TRUE)
}

.biohelper_stage_validation_summary <- function(value) {
  if (is.data.frame(value)) {
    list(rows = nrow(value), columns = names(value), types = vapply(value, typeof, character(1)))
  } else list(hits = value$hit_count, queries = value$query_count)
}

.biohelper_stage <- function(name, label, output, signature, run, validate,
                              checkpoint_dir, resume = TRUE, overwrite = FALSE,
                              adopt = TRUE) {
  metadata_file <- file.path(checkpoint_dir, paste0(name, ".rds"))
  metadata <- if (file.exists(metadata_file)) tryCatch(readRDS(metadata_file), error = function(e) NULL) else NULL
  key <- .biohelper_object_hash(signature)
  # A failed replacement must not discard the checkpoint of the completed file
  # still at its original path. This also allows returning to its parameters.
  previous <- if (identical(metadata$status, "complete")) metadata else metadata$previous_complete
  if (resume && !identical(metadata$status, "complete") && !is.null(previous) &&
      identical(previous$key, key) && identical(previous$output, .biohelper_file_fingerprint(output))) {
    metadata <- previous
    .biohelper_atomic_rds(metadata, metadata_file)
  }
  check <- function(path) tryCatch(list(value = validate(path), error = NULL), error = function(e) list(value = NULL, error = conditionMessage(e)))
  if (resume && !is.null(metadata) && identical(metadata$status, "complete") && identical(metadata$key, key) && file.exists(output)) {
    fingerprint <- .biohelper_file_fingerprint(output)
    if (identical(fingerprint, metadata$output)) {
      checked <- check(output)
      if (is.null(checked$error)) {
        message(label, ": SKIPPED (valid checkpoint)")
        return(list(value = checked$value, metadata = metadata, status = "skipped"))
      }
    }
  }
  # Older files are adopted only after full structural validation. Provenance
  # cannot be reconstructed from their contents and is explicitly marked unknown.
  if (resume && adopt && is.null(metadata) && file.exists(output) && file.info(output)$size > 0) {
    checked <- check(output)
    if (is.null(checked$error)) {
      warning(label, ": adopting existing file without checkpoint metadata; original parameters and database/software versions are unverified.", call. = FALSE)
      metadata <- list(schema = 1L, stage = name, status = "complete", key = key,
                       signature = signature, adopted = TRUE, provenance_verified = FALSE,
                       generation = .biohelper_object_hash(list(key, as.numeric(Sys.time()))),
                       completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       output = .biohelper_file_fingerprint(output),
                       validation = .biohelper_stage_validation_summary(checked$value))
      .biohelper_atomic_rds(metadata, metadata_file)
      message(label, ": SKIPPED (validated legacy file; provenance unverified)")
      return(list(value = checked$value, metadata = metadata, status = "adopted"))
    }
  }
  if (!resume && !overwrite && (file.exists(output) || file.exists(metadata_file))) {
    stop(label, ": existing results require resume=TRUE, a new output directory, or overwrite=TRUE.", call. = FALSE)
  }
  message(label, ": RUNNING")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(output), "."), tmpdir = dirname(output))
  on.exit(unlink(temporary))
  running <- list(schema = 1L, stage = name, status = "running", key = key, signature = signature,
                  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  if (!is.null(previous)) {
    previous$previous_complete <- NULL
    running$previous_complete <- previous
  }
  .biohelper_atomic_rds(running, metadata_file)
  analysis <- run(temporary)
  checked <- check(temporary)
  if (!is.null(checked$error)) stop(label, ": output validation failed: ", checked$error, call. = FALSE)
  .biohelper_atomic_publish(temporary, output, overwrite = overwrite, provenance = previous)
  complete <- c(running, list(completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                              generation = .biohelper_object_hash(list(key, temporary, as.numeric(Sys.time()))),
                              adopted = FALSE, provenance_verified = TRUE,
                              output = .biohelper_file_fingerprint(output),
                              validation = .biohelper_stage_validation_summary(checked$value),
                              analysis = if (is.list(analysis) && !is.data.frame(analysis)) analysis else NULL))
  complete$status <- "complete"
  complete$previous_complete <- NULL
  .biohelper_atomic_rds(complete, metadata_file)
  message(label, ": COMPLETED")
  list(value = checked$value, metadata = complete, status = "completed")
}

.biohelper_read_assignment_csv <- function(path, processed = FALSE, ids = NULL) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Assignment table is missing or empty.", call. = FALSE)
  parsed <- suppressWarnings(readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
                                            na = c("NA", ""), show_col_types = FALSE))
  # as.data.frame drops readr's parsing-problem attribute; inspect it first.
  if (nrow(readr::problems(parsed))) stop("Malformed assignment CSV.", call. = FALSE)
  x <- as.data.frame(parsed)
  required <- if (processed) c("ASVs", "taxonomy", "pident", "qcovs", "evalue", "bitscore", "sseqid") else c("asv", tolower(.biohelper_taxonomy_output_ranks()), .biohelper_blast_stats_columns())
  if (!all(required %in% names(x))) stop("Assignment CSV has missing or outdated columns.", call. = FALSE)
  if (!processed && !identical(names(x), required)) stop("Final assignment CSV has unexpected columns or column order.", call. = FALSE)
  id_column <- if (processed) "ASVs" else "asv"
  if (anyNA(x[[id_column]]) || any(!nzchar(x[[id_column]])) || anyDuplicated(x[[id_column]])) stop("Assignment CSV has missing or duplicate ASV identifiers.", call. = FALSE)
  if (!is.null(ids) && any(!x[[id_column]] %in% ids)) stop("Assignment CSV contains ASVs absent from the query FASTA.", call. = FALSE)
  if (!processed && !is.null(ids) && !identical(sort(x$asv), sort(ids))) stop("Final assignment CSV does not contain every query ASV.", call. = FALSE)
  for (field in intersect(c("pident", "qcovs", "evalue", "bitscore"), names(x))) {
    value <- suppressWarnings(as.numeric(x[[field]]))
    if (any(!is.na(x[[field]]) & (!is.finite(value) | value < 0)) ||
        (field %in% c("pident", "qcovs") && any(value > 100, na.rm = TRUE))) stop("Invalid assignment statistics in ", field, call. = FALSE)
    x[[field]] <- value
  }
  if (processed) {
    widths <- lengths(strsplit(as.character(ifelse(is.na(x$taxonomy), "Unknown", x$taxonomy)), ";", fixed = TRUE))
    if (any(!widths %in% c(1L, 8L, 9L))) stop("Malformed processed taxonomy strings.", call. = FALSE)
    if (any(widths == 1L & !is.na(x$taxonomy) & !trimws(x$taxonomy) %in% c("Unknown", "NA", ""))) stop("Malformed processed taxonomy strings.", call. = FALSE)
  } else .biohelper_validate_blastn_taxo_assignment_output(x)
  x
}

.biohelper_pad_assignments <- function(out, ids) {
  out <- out[match(ids, out$asv), , drop = FALSE]
  out$asv <- ids
  rownames(out) <- NULL
  out
}

.biohelper_adopt_legacy_final <- function(source, output, ids, configuration, checkpoint_dir, overwrite) {
  original_fingerprint <- .biohelper_file_fingerprint(source)
  legacy <- tryCatch(withCallingHandlers(as.data.frame(data.table::fread(source,
    colClasses = "character", na.strings = c("NA", ""), check.names = FALSE, showProgress = FALSE)),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)), error = function(e) NULL)
  if (is.null(legacy) || nrow(legacy) == 0L) return(NULL)
  if (!any(tolower(names(legacy)) %in% c(tolower(.biohelper_taxonomy_output_ranks()), "superkingdom", "taxonomy"))) return(NULL)
  out <- tryCatch(.biohelper_format_blastn_taxo_assignment_output(legacy), error = function(e) NULL)
  if (is.null(out) || anyNA(out$asv) || anyDuplicated(out$asv) || any(!out$asv %in% ids)) return(NULL)
  out <- .biohelper_pad_assignments(out, ids)
  warning("Adopting a structurally valid final table without stage metadata. Original parameters, reference versions and unavailable hit statistics/provenance cannot be verified; missing statistics remain NA. Source file is preserved.", call. = FALSE)
  temporary <- tempfile(".legacy_final_", tmpdir = dirname(output))
  on.exit(unlink(temporary))
  readr::write_csv(out, temporary, na = "NA")
  .biohelper_atomic_publish(temporary, output, overwrite)
  metadata <- list(schema = 1L, stage = "consensus", status = "complete", adopted_final = TRUE,
    provenance_verified = FALSE, key = .biohelper_object_hash(configuration), signature = configuration,
    source = original_fingerprint, output = .biohelper_file_fingerprint(output))
  .biohelper_atomic_rds(metadata, file.path(checkpoint_dir, "consensus.rds"))
  .biohelper_atomic_rds(list(complete = TRUE, adopted_final = TRUE, input_provenance_verified = FALSE,
    versions = configuration$versions, reference_database = configuration$database,
    taxonomy_database = configuration$taxonomy, parameters = configuration$parameters,
    query_count = length(ids), stages = list(consensus = metadata),
    limitation = "Legacy parameters, software/reference versions and missing alignment statistics cannot be verified."),
    file.path(dirname(output), "blastn_taxo_assignment_provenance.rds"))
  message("Consensus merging: SKIPPED (adopted final; provenance unverified)")
  out
}

.biohelper_blast_versions <- function() {
  packages <- c("biohelper", "taxizedb", "taxonomizr", "DBI", "RSQLite", "readr", "data.table", "janitor")
  c(list(R = as.character(getRversion())), stats::setNames(lapply(packages, function(x) as.character(utils::packageVersion(x))), packages))
}

.biohelper_blast_db_identity <- function(db) {
  candidates <- c(db, Sys.glob(paste0(db, ".*")))
  if (!any(file.exists(candidates))) {
    directories <- strsplit(Sys.getenv("BLASTDB"), .Platform$path.sep, fixed = TRUE)[[1]]
    candidates <- unlist(lapply(directories[nzchar(directories)], function(directory) c(file.path(directory, db), Sys.glob(paste0(file.path(directory, db), ".*")))))
  }
  candidates <- as.character(candidates)
  candidates <- sort(unique(candidates[file.exists(candidates) & !dir.exists(candidates)]))
  list(reference = db, files = lapply(candidates, .biohelper_file_fingerprint, content = FALSE),
       identity_verified = length(candidates) > 0L,
       note = "Reference database identity uses file names, sizes and modification times, not full hashes of large sequence databases.")
}

.biohelper_query_ids <- function(path) {
  connection <- file(path, "r")
  on.exit(close(connection))
  ids <- character()
  repeat {
    lines <- readLines(connection, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    headers <- lines[startsWith(lines, ">")]
    ids <- c(ids, sub("[[:space:]].*$", "", substring(headers, 2L)))
  }
  if (!length(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("`queries` must be a FASTA file with non-empty, unique feature IDs.", call. = FALSE)
  sort(ids)
}

.biohelper_blast_executable <- function(path) {
  executable <- if (file.exists(path)) path else unname(Sys.which(path))
  if (!nzchar(executable) || !file.exists(executable)) stop("BLAST executable not found: ", path, call. = FALSE)
  result <- tryCatch(system2(executable, "-version", stdout = TRUE, stderr = TRUE), error = function(e) stop("Cannot run BLAST executable: ", conditionMessage(e), call. = FALSE))
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0L) stop("BLAST version check failed: ", paste(result, collapse = " "), call. = FALSE)
  list(path = normalizePath(executable), version = result[grepl("blastn", result)][1L])
}

.biohelper_validate_assignment_options <- function(options) {
  # The user-supplied option string is retained for API compatibility. Quote
  # controlled paths separately; disallow options that could replace outputs.
  if (!is.character(options) || length(options) != 1L || is.na(options) ||
      grepl("[\n\r;&|`<>]", options) || grepl("\\$\\(", options)) {
    stop("BLAST options must be a single option string without shell control characters.", call. = FALSE)
  }
  if (grepl("(^|[[:space:]])-(out|outfmt|query|db|task|num_threads)([[:space:]=]|$)", options)) {
    stop("Set BLAST inputs, database, task, threads and output through the function arguments, not the option string.", call. = FALSE)
  }
  invisible(options)
}

.biohelper_run_assignment_search <- function(executable, queries, db, task, options, nthreads, output, log) {
  .biohelper_validate_assignment_options(options)
  format <- "6 qseqid qlen pident sseqid sgi evalue bitscore staxids sscinames sblastnames qcovs qcovhsp"
  temporary_log <- tempfile("blast_log_", tmpdir = dirname(output))
  on.exit(unlink(temporary_log))
  arguments <- c("-db", shQuote(db), "-query", shQuote(queries), "-num_threads", nthreads,
                 options, "-task", task, "-out", shQuote(output), "-outfmt", shQuote(format))
  status <- tryCatch(system2(executable, arguments, stdout = temporary_log, stderr = temporary_log),
                     error = function(e) { writeLines(conditionMessage(e), temporary_log); 127L })
  text <- readLines(temporary_log, warn = FALSE)
  cat(paste0("\n", format(Sys.time(), tz = "UTC", usetz = TRUE), " ", task, " exit=", status, "\n"),
      paste(text, collapse = "\n"), "\n", file = log, append = TRUE)
  if (status != 0L) stop(task, " failed (exit ", status, "). See ", log, ". ", paste(utils::tail(text, 5L), collapse = " "), call. = FALSE)
  if (!file.exists(output)) stop(task, " returned success but did not create an output file. See ", log, call. = FALSE)
}
