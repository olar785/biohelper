.biohelper_blast_http_state <- new.env(parent = emptyenv())
.biohelper_blast_http_state$last_request <- -Inf
.biohelper_blast_now <- function() as.numeric(Sys.time())
.biohelper_blast_wait <- function(seconds) Sys.sleep(seconds)

.biohelper_blast_http <- function(parameters) {
  if (!requireNamespace("curl", quietly = TRUE)) stop("Install the suggested package `curl` for remote BLAST.", call. = FALSE)
  elapsed <- .biohelper_blast_now() - .biohelper_blast_http_state$last_request
  if (elapsed < 10) .biohelper_blast_wait(10 - elapsed)
  encoded <- paste(vapply(names(parameters), function(name) {
    paste0(utils::URLencode(name, reserved = TRUE), "=",
           utils::URLencode(as.character(parameters[[name]]), reserved = TRUE))
  }, character(1)), collapse = "&")
  handle <- curl::new_handle()
  curl::handle_setopt(handle, postfields = encoded, connecttimeout = 30, timeout = 120,
                      useragent = "biohelper NCBI BLAST checker")
  curl::handle_setheaders(handle, "Content-Type" = "application/x-www-form-urlencoded")
  .biohelper_blast_http_state$last_request <- .biohelper_blast_now()
  result <- curl::curl_fetch_memory("https://blast.ncbi.nlm.nih.gov/Blast.cgi", handle)
  if (result$status_code != 200L) stop("NCBI BLAST HTTP status ", result$status_code, call. = FALSE)
  rawToChar(result$content)
}

.biohelper_parse_blast_json <- function(text, query_ids) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install the suggested package `jsonlite` for remote BLAST.", call. = FALSE)
  document <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  reports <- document$BlastOutput2
  if (is.null(reports) && !is.null(document$report)) reports <- list(document)
  if (is.null(reports)) stop("NCBI response does not contain BlastOutput2 reports.", call. = FALSE)
  rows <- list()
  seen <- character()
  for (entry in reports) {
    if (!is.null(entry$error)) stop("NCBI BLAST error: ", paste(unlist(entry$error), collapse = "; "), call. = FALSE)
    search <- entry$report$results$search
    if (is.null(search)) stop("NCBI BLAST report is missing query results.", call. = FALSE)
    query <- if (!is.null(search$query_title)) sub("[[:space:]].*$", "", search$query_title) else search$query_id
    if (!query %in% query_ids || query %in% seen) stop("Unexpected or duplicate query identifier in NCBI BLAST response: ", query, call. = FALSE)
    seen <- c(seen, query)
    query_length <- as.numeric(search$query_len)
    if (length(query_length) != 1L || !is.finite(query_length) || query_length <= 0) stop("Invalid NCBI query length.", call. = FALSE)
    for (hit in search$hits) {
      descriptions <- hit$description
      if (!length(descriptions) || !length(hit$hsps)) next
      subject <- descriptions[[1]]$accession
      if (is.null(subject)) subject <- descriptions[[1]]$id
      taxa <- unique(vapply(descriptions, function(x) if (is.null(x$taxid)) "N/A" else as.character(x$taxid), character(1)))
      scientific <- unique(vapply(descriptions, function(x) if (is.null(x$sciname)) "N/A" else x$sciname, character(1)))
      # qcovs is the union of query intervals; screening uses each HSP separately.
      intervals <- t(vapply(hit$hsps, function(hsp) sort(c(hsp$query_from, hsp$query_to)), numeric(2)))
      intervals <- intervals[order(intervals[, 1]), , drop = FALSE]
      total <- 0
      right <- 0
      for (i in seq_len(nrow(intervals))) {
        total <- total + max(0, intervals[i, 2] - max(right, intervals[i, 1] - 1))
        right <- max(right, intervals[i, 2])
      }
      for (hsp in hit$hsps) {
        rows[[length(rows) + 1L]] <- data.frame(
          qseqid = query, qlen = query_length, pident = 100 * hsp$identity / hsp$align_len,
          sseqid = subject, sgi = 0, evalue = hsp$evalue, bitscore = hsp$bit_score,
          staxids = paste(taxa, collapse = ";"), sscinames = paste(scientific, collapse = ";"),
          sblastnames = "N/A", qcovs = 100 * total / query_length,
          qcovhsp = 100 * (abs(hsp$query_to - hsp$query_from) + 1) / query_length,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!setequal(seen, query_ids)) stop("NCBI BLAST response is missing query results; no partial result was cached.", call. = FALSE)
  if (length(rows)) .biohelper_validate_blast_hits(do.call(rbind, rows)) else .biohelper_empty_blast_hits()
}

.biohelper_run_remote_blast <- function(query, output, log, control) {
  if (!requireNamespace("curl", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Remote BLAST requires the suggested packages `curl` and `jsonlite`.", call. = FALSE)
  }
  metadata_file <- paste0(output, ".request.rds")
  query_text <- paste(readLines(query, warn = FALSE), collapse = "\n")
  query_ids <- sub("[[:space:]].*$", "", sub("^>", "", grep("^>", readLines(query, warn = FALSE), value = TRUE)))
  if (!control$refresh && file.exists(metadata_file)) {
    metadata <- readRDS(metadata_file)
    if (!identical(metadata$query, query_text)) stop("BLAST request cache does not match the submitted queries.", call. = FALSE)
    .biohelper_blast_message(control, "Resuming NCBI BLAST request ", metadata$rid, ".")
  } else {
    parameters <- list(CMD = "Put", PROGRAM = "blastn", DATABASE = control$db,
                       QUERY = query_text, EXPECT = control$evalue, HITLIST_SIZE = control$max_hits, tool = "biohelper")
    if (control$task == "megablast") parameters$MEGABLAST <- "on"
    if (!is.null(control$email)) parameters$email <- control$email
    response <- .biohelper_blast_http(parameters)
    writeLines(response, log)
    match <- regexec("RID[[:space:]]*=[[:space:]]*([A-Za-z0-9-]+)", response)
    rid <- regmatches(response, match)[[1]]
    if (length(rid) < 2L) stop("NCBI did not return a BLAST request ID. See ", log, call. = FALSE)
    metadata <- list(rid = rid[[2]], query = query_text, last_poll = .biohelper_blast_now())
    saveRDS(metadata, metadata_file)
    .biohelper_blast_message(control, "NCBI BLAST request ", metadata$rid, " submitted; request metadata saved for resume.")
  }
  start <- .biohelper_blast_now()
  repeat {
    delay <- max(0, control$poll_seconds - (.biohelper_blast_now() - metadata$last_poll))
    if (.biohelper_blast_now() - start + delay >= control$timeout_seconds) {
      stop("NCBI BLAST request ", metadata$rid, " is still pending. Rerun with the same cache directory to resume.", call. = FALSE)
    }
    if (delay > 0) .biohelper_blast_wait(delay)
    metadata$last_poll <- .biohelper_blast_now()
    saveRDS(metadata, metadata_file)
    response <- .biohelper_blast_http(list(CMD = "Get", RID = metadata$rid, FORMAT_TYPE = "JSON2_S", tool = "biohelper"))
    writeLines(response, log)
    if (grepl("Status[[:space:]]*=[[:space:]]*WAITING", response)) {
      .biohelper_blast_message(control, "NCBI BLAST request ", metadata$rid,
        " is still waiting (", round(.biohelper_blast_now() - start, 1), " s elapsed this call).")
      next
    }
    if (grepl("Status[[:space:]]*=[[:space:]]*(FAILED|UNKNOWN)", response)) {
      writeLines(response, log)
      stop("NCBI BLAST request ", metadata$rid, " failed or expired. See ", log, call. = FALSE)
    }
    writeLines(response, paste0(output, ".json"))
    hits <- .biohelper_parse_blast_json(response, query_ids)
    if (nrow(hits)) utils::write.table(hits, output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "N/A") else file.create(output)
    return(invisible(output))
  }
}
