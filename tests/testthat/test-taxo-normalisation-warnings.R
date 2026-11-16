capture_taxonomy_warnings <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

test_that("lineage-resolved homonyms produce no stale lookup warning in either entry point", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = c("a", "b"), broad = "Metazoa", fine = "Prunella")
  result <- capture_taxonomy_warnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check",
    blast_options = list(verbose = FALSE)
  ))
  expect_length(result$warnings, 0L)
  expect_equal(result$value$family, rep("Prunellidae", 2))
  expect_equal(nrow(attr(result$value, "blast_check")$diagnostics), 0L)
  merged <- capture_taxonomy_warnings(taxo_merge(
    list(silva = input, pr2 = input), database, addExtra = FALSE,
    ambiguity = "blast_check", blast_options = list(verbose = FALSE)
  ))
  expect_length(merged$warnings, 0L)
  expect_equal(merged$value$family, rep("Prunellidae", 2))
})

test_that("unresolved homonyms retain structured ambiguity warnings for every policy", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "unresolved", broad = "Eukaryota", fine = "Prunella")
  for (policy in c("first", "lca", "unassigned")) {
    result <- capture_taxonomy_warnings(taxo_normalisation(input, database,
      addExtra = FALSE, ambiguity = policy))
    expect_length(result$warnings, 1L)
    expect_s3_class(result$warnings[[1]], "biohelper_taxonomy_ambiguity")
    expect_identical(result$warnings[[1]]$feature_ids, "unresolved")
    expect_false(grepl("Collapsing with commas", conditionMessage(result$warnings[[1]]), fixed = TRUE))
  }
  # Without any distinguishing broader clue, the first policy leaves NA.
  result <- capture_taxonomy_warnings(taxo_normalisation(input[, c("feature_id", "fine")],
    database, addExtra = FALSE))
  expect_length(result$warnings, 1L)
  expect_s3_class(result$warnings[[1]], "biohelper_taxonomy_ambiguity")
  expect_match(conditionMessage(result$warnings[[1]]), "Left 1 feature(s) unassigned", fixed = TRUE)
  expect_true(is.na(result$value$genus))
})

test_that("BLAST-resolved names emit no lookup warning and insufficient evidence still warns", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", broad = "Eukaryota", fine = "Prunella")
  hits <- data.frame(qseqid = "asv1", sseqid = c("a", "b"), staxids = "175132",
    pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
  result <- capture_taxonomy_warnings(taxo_normalisation(input, database, addExtra = FALSE,
    ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(min_taxa = 1L, verbose = FALSE)))
  expect_length(result$warnings, 0L)
  expect_equal(result$value$family, "Prunellidae")
  result <- capture_taxonomy_warnings(taxo_normalisation(input, database, addExtra = FALSE,
    ambiguity = "blast_check", blast_hits = hits, blast_options = list(verbose = FALSE)))
  expect_equal(length(result$warnings), 2L)
  expect_true(any(vapply(result$warnings, inherits, logical(1), "biohelper_blast_review")))
  expect_true(any(vapply(result$warnings, inherits, logical(1), "biohelper_taxonomy_ambiguity")))
})

test_that("unrelated lookup warnings and errors are not suppressed", {
  database <- ncbi_prunella_database()
  original <- taxonomizr::getId
  testthat::local_mocked_bindings(getId = function(...) {
    warning(structure(list(message = "Other lookup issue", call = NULL),
      class = c("test_lookup_warning", "warning", "condition")))
    original(...)
  }, .package = "taxonomizr")
  input <- data.frame(feature_id = "asv1", broad = "Metazoa", fine = "Prunella")
  result <- capture_taxonomy_warnings(taxo_normalisation(input, database, addExtra = FALSE))
  expect_length(result$warnings, 1L)
  expect_s3_class(result$warnings[[1]], "test_lookup_warning")
  expect_equal(conditionMessage(result$warnings[[1]]), "Other lookup issue")
  expect_equal(result$value$family, "Prunellidae")
})

test_that("lookup failures still stop normalisation", {
  testthat::local_mocked_bindings(getId = function(...) stop("Lookup failed"), .package = "taxonomizr")
  input <- data.frame(feature_id = "asv1", fine = "Prunella")
  expect_error(taxo_normalisation(input, "unused", addExtra = FALSE), "Lookup failed")
})
