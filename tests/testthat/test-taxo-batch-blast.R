batch_taxonomy_input <- function(ids) {
  data.frame(feature_id = ids, broad = "Eukaryota", fine = "Prunella")
}

test_that("progress reports distinct ambiguous ASVs rather than summing duplicate table entries", {
  database <- ncbi_prunella_database()
  messages <- character()
  withCallingHandlers(suppressWarnings(taxo_merge(
    list(batch_taxonomy_input(c("a", "b")), batch_taxonomy_input("a")), database,
    addExtra = FALSE, ambiguity = "blast_check", blast_options = list(fallback = "lca")
  )), message = function(m) {
    messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  expect_true(any(grepl("Ambiguous ASVs found: 2 distinct across all input tables", messages, fixed = TRUE)))
  expect_true(any(grepl("3 selected assignments across 2 tables; 2 distinct ASVs", messages, fixed = TRUE)))
  expect_true(any(grepl("2 ambiguous ASVs", messages, fixed = TRUE)))
  expect_true(any(grepl("2 selected ASVs have no sequence", messages, fixed = TRUE)))
})

test_that("one pooled search covers disjoint tables and deduplicates DNA without pooling taxonomy evidence", {
  database <- ncbi_prunella_database()
  tables <- list(silva = batch_taxonomy_input(c("bird", "clone")),
    pr2 = batch_taxonomy_input(c("plant", "bird2")),
    other = batch_taxonomy_input(c("plant_clone", "absent")))
  original <- tables
  dna <- c(bird = "ACGTACGT", clone = "ACGTACGT", plant = "TTTTTTTT",
    bird2 = "GGGGGGGG", plant_clone = "TTTTTTTT")
  cache <- tempfile("pooled-blast-")
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  searches <- 0L
  query_count <- 0L
  testthat::local_mocked_bindings(.biohelper_run_blast = function(query, output, log, control) {
    searches <<- searches + 1L
    queries <- Biostrings::readDNAStringSet(query)
    query_count <<- length(queries)
    records <- lapply(seq_along(queries), function(i) {
      data.frame(qseqid = names(queries)[[i]], qlen = 8, pident = 99,
        sseqid = c("reference_a", "reference_b"), sgi = 0, evalue = 1e-30,
        bitscore = 100, staxids = if (as.character(queries[[i]]) == "TTTTTTTT") "39176" else "175132",
        sscinames = "Prunella", sblastnames = "test", qcovs = 100, qcovhsp = 100)
    })
    utils::write.table(do.call(rbind, records), output, quote = FALSE, row.names = FALSE,
      col.names = FALSE, sep = "\t")
  })
  args <- list(df_list = tables, sqlFile = database, sequences = dna, addExtra = FALSE,
    ambiguity = "blast_check", db = "/database/nt", blast_options = list(
      fallback = "lca", min_taxa = 1L, cache_dir = cache, verbose = FALSE))
  result <- suppressWarnings(do.call(taxo_merge, args))
  expect_identical(tables, original)
  expect_equal(searches, 1L)
  expect_equal(query_count, 3L)
  expect_equal(result$family[match(names(dna), result$feature_id)],
    c("Prunellidae", "Prunellidae", "Lamiaceae", "Prunellidae", "Lamiaceae"))
  expect_true(is.na(result$genus[result$feature_id == "absent"]))
  reports <- attr(result, "normalisation")
  expect_named(reports, names(tables))
  expect_equal(reports$other$blast$diagnostics$status, c("blast_resolved", "missing_sequence"))
  expect_equal(reports$pr2$blast$hit_source, "blastn")
  expect_equal(reports$pr2$blast$diagnostics$feature_id, c("plant", "bird2"))
  again <- suppressWarnings(do.call(taxo_merge, args))
  expect_identical(again, result)
  expect_equal(searches, 1L)
  args$blast_options$refresh <- TRUE
  suppressWarnings(do.call(taxo_merge, args))
  expect_equal(searches, 2L) # One fresh batch, not one refresh per table.
})

test_that("matching repeated taxonomy profiles does not merge ASV-specific sequence evidence", {
  database <- ncbi_prunella_database()
  input <- batch_taxonomy_input(c("a", "b", "c", "d"))
  matches <- character()
  original_get_id <- taxonomizr::getId
  testthat::local_mocked_bindings(getId = function(taxa, ...) {
    matches <<- c(matches, taxa)
    original_get_id(taxa, ...)
  }, .package = "taxonomizr")
  result <- suppressWarnings(taxo_normalisation(input, database, addExtra = FALSE,
    ambiguity = "lca"))
  expect_equal(result$feature_id, input$feature_id)
  expect_equal(sum(matches == "Prunella"), 1L)
  expect_equal(result$domain, rep("Eukaryota", 4))
  expect_equal(nrow(attr(result, "ambiguity")), 4L)
})

test_that("a global feature subset can span tables with different ASV sets", {
  database <- ncbi_prunella_database()
  searches <- 0L
  testthat::local_mocked_bindings(.biohelper_blast_search = function(sequences, control) {
    searches <<- searches + 1L
    expect_setequal(names(sequences), c("a", "b"))
    biohelper:::.biohelper_empty_blast_hits()
  })
  tables <- list(batch_taxonomy_input("a"), batch_taxonomy_input("b"))
  result <- suppressWarnings(taxo_merge(tables, database, addExtra = FALSE,
    ambiguity = "blast_check", sequences = c(a = "ACGT", b = "GGGG"), db = "nt",
    blast_options = list(feature_ids = c("a", "b"), fallback = "lca", verbose = FALSE)))
  expect_equal(searches, 1L)
  expect_true(all(is.na(result$genus)))
  expect_error(suppressWarnings(taxo_merge(tables, database, ambiguity = "blast_check",
    blast_options = list(feature_ids = "unknown", verbose = FALSE))), "Unknown BLAST feature IDs")
})

test_that("empty evidence and no sequence files do not start a local process", {
  database <- ncbi_prunella_database()
  testthat::local_mocked_bindings(.biohelper_run_blast = function(...) stop("Unexpected search"))
  result <- suppressWarnings(taxo_merge(list(batch_taxonomy_input("a")), database,
    ambiguity = "blast_check", addExtra = FALSE,
    blast_options = list(fallback = "lca", verbose = FALSE)))
  expect_equal(attr(result, "normalisation")[[1]]$blast$diagnostics$status, "missing_sequence")
  unique_name <- data.frame(feature_id = "b", name = "Lamiaceae")
  result <- taxo_merge(list(unique_name), database, addExtra = FALSE, ambiguity = "blast_check",
    sequences = "/no/FASTA/needed", blast_options = list(verbose = FALSE))
  expect_equal(result$family, "Lamiaceae")
  expect_equal(nrow(attr(result, "normalisation")[[1]]$blast$diagnostics), 0L)
})

test_that("the pooled remote query limit is checked before any submission", {
  database <- ncbi_prunella_database()
  testthat::local_mocked_bindings(.biohelper_run_blast = function(...) stop("Unexpected remote submission"))
  expect_error(suppressWarnings(taxo_merge(
    list(batch_taxonomy_input("a"), batch_taxonomy_input("b")), database,
    ambiguity = "blast_check", addExtra = FALSE, sequences = c(a = "ACGT", b = "GGGG"),
    blast_options = list(remote = TRUE, max_remote_queries = 1L, cache_dir = tempfile(), verbose = FALSE)
  )), "2 uncached remote queries")
})

test_that("audit mode keeps per-input candidates when sequence evidence is shared", {
  database <- ncbi_prunella_database()
  searches <- 0L
  testthat::local_mocked_bindings(.biohelper_blast_search = function(sequences, control) {
    searches <<- searches + 1L
    do.call(rbind, lapply(names(sequences), function(id) data.frame(
      qseqid = id, sseqid = c("a", "b"), staxids = "175132", pident = 99,
      qcovhsp = 100, evalue = 1e-30, bitscore = 100)))
  })
  tables <- list(plant = data.frame(feature_id = "plant", name = "Lamiaceae"),
    bird = data.frame(feature_id = "bird", name = "Prunellidae"))
  result <- suppressWarnings(taxo_merge(tables, database, addExtra = FALSE,
    ambiguity = "blast_check", sequences = c(plant = "ACGT", bird = "ACGT"), db = "nt",
    blast_options = list(scope = "all", min_taxa = 1L, conflict_action = "unassign", verbose = FALSE)))
  expect_equal(searches, 1L)
  expect_equal(result$feature_id, "bird")
  reports <- attr(result, "normalisation")
  expect_equal(reports$plant$blast$diagnostics$status, "blast_conflict")
  expect_equal(reports$bird$blast$diagnostics$status, "blast_supported")
})

test_that("progress settings are validated before matching or searching", {
  for (control in list(list(verbose = NA), list(verbose = 1), list(progress_seconds = 0),
    list(progress_seconds = NA), list(progress_seconds = "30"))) {
    expect_error(biohelper:::.biohelper_blast_options(control), "Invalid BLAST")
  }
  defaults <- biohelper:::.biohelper_blast_options(list())
  expect_true(defaults$verbose)
  expect_equal(defaults$progress_seconds, 30)
})

test_that("local searches show elapsed-time heartbeats and preserve paths with spaces", {
  skip_on_os("windows")
  directory <- tempfile("heartbeat with spaces ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  executable <- file.path(directory, "fake blastn")
  writeLines(c("#!/bin/sh", "sleep 0.15", 'while [ "$#" -gt 0 ]; do',
    '  if [ "$1" = "-out" ]; then shift; : > "$1"; fi', "  shift", "done"), executable)
  Sys.chmod(executable, "0755")
  control <- biohelper:::.biohelper_blast_options(list(progress_seconds = 0.02))
  control$db <- file.path(directory, "database prefix")
  control$blastapp_path <- executable
  output <- file.path(directory, "output hits.tsv")
  messages <- character()
  withCallingHandlers(biohelper:::.biohelper_run_blast("query path", output,
    file.path(directory, "search.log"), control), message = function(m) {
      messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  expect_true(file.exists(output))
  expect_true(any(grepl("BLAST still running", messages, fixed = TRUE)))
  control$verbose <- FALSE
  expect_message(biohelper:::.biohelper_run_blast("query path", output,
    file.path(directory, "quiet.log"), control), NA)
  writeLines(c("#!/bin/sh", "exit 7"), executable)
  expect_error(biohelper:::.biohelper_run_blast("query path", output,
    file.path(directory, "failure.log"), control), "BLAST failed")
})

test_that("an interrupted local search cleans up its child process", {
  killed <- FALSE
  fake <- list(is_alive = function() TRUE, wait = function(...) stop("interrupted"),
    kill = function() { killed <<- TRUE })
  testthat::local_mocked_bindings(.biohelper_start_blast_process = function(...) fake)
  control <- biohelper:::.biohelper_blast_options(list(verbose = FALSE))
  control$db <- "nt"
  control$blastapp_path <- tempfile()
  file.create(control$blastapp_path)
  on.exit(unlink(control$blastapp_path), add = TRUE)
  expect_error(biohelper:::.biohelper_run_blast("query", "output", "log", control), "interrupted")
  expect_true(killed)
})
