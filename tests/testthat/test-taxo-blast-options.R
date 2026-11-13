test_that("the public BLAST interface defaults to local and matches existing argument names", {
  defaults <- biohelper:::.biohelper_blast_options(list())
  expect_false(defaults$remote)
  expect_equal(defaults$scope, "ambiguous")
  expect_equal(defaults$fallback, "first")
  expect_equal(defaults$nthreads, 1L)
  expect_true(all(c("db", "blastapp_path", "blast_options") %in% names(formals(taxo_normalisation))))
  expect_false("blast_control" %in% names(formals(taxo_normalisation)))
  expect_true(all(c("db", "blastapp_path") %in% names(formals(blastn_taxo_assignment))))

  control <- defaults
  control$db <- NULL
  expect_error(biohelper:::.biohelper_run_blast("query", "out", "log", control), "Supply `db`")
  control$db <- "/databases/nt"
  control$blastapp_path <- "/nonexistent/blastn"
  expect_error(biohelper:::.biohelper_run_blast("query", "out", "log", control), "blastapp_path")
})

test_that("local database and executable options reach the search backend", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", a = "Eukaryota", b = "Prunella")
  observed <- NULL
  testthat::local_mocked_bindings(.biohelper_blast_search = function(sequences, control) {
    observed <<- control
    expect_identical(sequences, c(asv1 = "ACGT"))
    data.frame(qseqid = "asv1", sseqid = c("a", "b"), staxids = "175132",
               pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
  })
  output <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", sequences = c(asv1 = "ACGT"),
    db = "/database with spaces/nt", blastapp_path = "/tools/blastn",
    blast_options = list(fallback = "lca", min_taxa = 1L, nthreads = 3L), return_details = TRUE
  ))
  expect_false(observed$remote)
  expect_equal(observed$db, "/database with spaces/nt")
  expect_equal(observed$blastapp_path, "/tools/blastn")
  expect_equal(observed$nthreads, 3L)
  expect_equal(output$taxonomy$family, "Prunellidae")
  expect_equal(output$blast$diagnostics$status, "blast_resolved")

  output <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", sequences = c(asv1 = "ACGT"),
    blast_options = list(remote = TRUE, min_taxa = 1L), return_details = TRUE
  ))
  expect_true(observed$remote)
  expect_equal(observed$db, "core_nt")
  expect_equal(output$taxonomy$family, "Prunellidae")
})

test_that("missing sequence IDs automatically use LCA without manual intersections", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = c("present", "absent"), a = "Eukaryota", b = "Prunella")
  original <- input
  searches <- 0L
  testthat::local_mocked_bindings(.biohelper_blast_search = function(sequences, control) {
    searches <<- searches + 1L
    expect_equal(names(sequences), "present")
    data.frame(qseqid = "present", sseqid = c("a", "b"), staxids = "175132",
               pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
  })
  fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">present description", "ACGT", ">unrelated", "AAAA"), fasta)
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", sequences = fasta,
    db = "/databases/nt", blast_options = list(fallback = "lca", min_taxa = 1L), return_details = TRUE
  ))
  expect_equal(searches, 1L)
  expect_identical(input, original)
  expect_equal(result$blast$diagnostics$status, c("blast_resolved", "missing_sequence"))
  expect_equal(result$taxonomy$family[[1]], "Prunellidae")
  expect_true(is.na(result$taxonomy$family[[2]]))
  expect_equal(result$ambiguity$lca_id, "2759")
  expect_equal(result$blast$diagnostics$final_id, c("175132", "2759"))

  # No actual search means no local database or executable is needed.
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check",
    blast_options = list(fallback = "lca"), return_details = TRUE
  ))
  expect_equal(searches, 1L)
  expect_true(all(result$blast$diagnostics$status == "missing_sequence"))
  expect_equal(result$ambiguity$lca_id, c("2759", "2759"))
})
