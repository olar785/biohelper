test_that("one ambiguity selector offers first, LCA and unassigned policies", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", phylum = "Eukaryota", genus = "Prunella")
  default <- suppressWarnings(taxo_normalisation(input, database, addExtra = FALSE))
  first <- suppressWarnings(taxo_normalisation(input, database, addExtra = FALSE, ambiguity = "first"))
  expect_identical(default, first)
  expect_null(attr(first, "ambiguity"))
  expect_equal(first$family, "Lamiaceae")
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "lca", return_details = TRUE
  ))
  expect_equal(result$taxonomy$domain, "Eukaryota")
  expect_true(is.na(result$taxonomy$family))
  expect_true(is.na(result$taxonomy$genus))
  expect_equal(result$ambiguity$lca_id, "2759")
  expect_equal(result$ambiguity$lca_name, "Eukaryota")
  expect_equal(result$ambiguity$lca_rank, "domain")
  expect_null(result$blast)
  unassigned <- suppressWarnings(taxo_normalisation(input, database, addExtra = FALSE, ambiguity = "unassigned"))
  expect_equal(unassigned$domain, "Unknown")
  expect_true(is.na(unassigned$genus))
  expect_equal(attr(unassigned, "ambiguity")$action, "unassigned")
  expect_false("blast_check" %in% names(formals(taxo_normalisation)))
})

test_that("LCA uses ancestry IDs for nested candidates and preserves omitted rank information", {
  database <- ncbi_prunella_database()
  family <- taxonomizr::getId("Lamiaceae", database)
  result <- biohelper:::.biohelper_resolve_ambiguity(
    data.frame(feature_id = "asv1", id = "39176"), list(c("39176", family)), 1L, "lca", database
  )
  expect_equal(result$ids, family)
  expect_equal(result$diagnostics$lca_name, "Lamiaceae")
  expect_equal(result$diagnostics$lca_rank, "family")
  input <- data.frame(feature_id = "asv1", phylum = NA_character_, genus = "Prunella")
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "lca", ranks = c("family", "genus"), return_details = TRUE
  ))
  expect_named(result$taxonomy, c("feature_id", "family", "genus"))
  expect_true(all(is.na(result$taxonomy[, c("family", "genus")])))
  expect_equal(result$ambiguity$lca_name, "Eukaryota")
})

test_that("BLAST is an ambiguity method with configurable fallback", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", phylum = "Eukaryota", genus = "Prunella")
  hits <- data.frame(qseqid = "asv1", sseqid = c("record1", "record2"), staxids = "175132",
                     pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
  resolved <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(min_taxa = 1L, fallback = "lca"), return_details = TRUE
  ))
  expect_equal(resolved$taxonomy$family, "Prunellidae")
  expect_equal(resolved$blast$diagnostics$status, "blast_resolved")
  expect_equal(nrow(resolved$ambiguity), 0L)
  weak <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(fallback = "lca"), return_details = TRUE
  ))
  expect_equal(weak$blast$diagnostics$status, "insufficient_records")
  expect_equal(weak$taxonomy$domain, "Eukaryota")
  expect_true(is.na(weak$taxonomy$genus))
  expect_equal(weak$ambiguity$action, "lca")
  expect_equal(weak$ambiguity$method, "blast_check")
  expect_equal(weak$blast$diagnostics$final_id, "2759")
  unassigned <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(fallback = "unassigned"), return_details = TRUE
  ))
  expect_equal(unassigned$taxonomy$domain, "Unknown")
  expect_true(is.na(unassigned$blast$diagnostics$final_id))
  expect_error(biohelper:::.biohelper_blast_options(list(fallback = "blast_check")), "fallback")
})

test_that("LCA fallback cannot reinstate an explicit BLAST conflict rejection", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", phylum = "Eukaryota", genus = "Prunella")
  testthat::local_mocked_bindings(.biohelper_check_blast = function(matches, candidates, sqlFile, mode, sequences, hits, control, ...) {
    list(ids = NA_character_, diagnostics = data.frame(feature_id = "asv1", status = "blast_conflict", final_id = NA_character_),
         hits = data.frame(), settings = control)
  })
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check",
    blast_options = list(fallback = "lca", conflict_action = "unassign"), return_details = TRUE
  ))
  expect_equal(result$taxonomy$domain, "Unknown")
  expect_true(is.na(result$taxonomy$genus))
  expect_equal(nrow(result$ambiguity), 0L)
})
