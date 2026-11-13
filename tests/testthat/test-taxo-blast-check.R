blast_test_paths <- function() {
  ancestors <- list(`1` = "1", `10` = c("10", "1"), `11` = c("11", "10", "1"),
                    `20` = c("20", "1"), `101` = c("101", "10", "1"),
                    `102` = c("102", "10", "1"), `111` = c("111", "11", "10", "1"),
                    `112` = c("112", "11", "10", "1"), `201` = c("201", "20", "1"),
                    `202` = c("202", "20", "1"), `301` = c("301", "1"), `302` = c("302", "1"))
  lapply(ancestors, function(ids) list(ids = ids, complete = TRUE))
}

blast_test_hits <- function(taxids = c("101", "102")) {
  data.frame(qseqid = "asv1", sseqid = paste0("record", seq_along(taxids)), staxids = taxids,
             pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
}

test_that("sequence evidence resolves a branch or a supported nested candidate", {
  paths <- blast_test_paths()
  testthat::local_mocked_bindings(.biohelper_ncbi_id_paths = function(ids, sqlFile) paths[ids])
  matches <- data.frame(feature_id = "asv1", id = "20")
  result <- biohelper:::.biohelper_check_blast(matches, list(c("10", "20")), "unused", "ambiguous",
                                             NULL, blast_test_hits(), list())
  expect_equal(result$ids, "10")
  expect_equal(result$diagnostics$status, "blast_resolved")
  expect_equal(result$diagnostics$n_taxa, 2L)
  result <- biohelper:::.biohelper_check_blast(matches, list(c("10", "11")), "unused", "ambiguous",
                                             NULL, blast_test_hits(c("111", "112")), list())
  expect_equal(result$ids, "11")
})

test_that("mixed, missing and weak evidence cannot silently choose a candidate", {
  paths <- blast_test_paths()
  testthat::local_mocked_bindings(.biohelper_ncbi_id_paths = function(ids, sqlFile) paths[ids])
  matches <- data.frame(feature_id = "asv1", id = "20")
  check <- function(hits, control = list()) biohelper:::.biohelper_check_blast(
    matches, list(c("10", "20")), "unused", "ambiguous", NULL, hits, control
  )
  result <- check(blast_test_hits(c("101", "201")))
  expect_equal(result$diagnostics$status, "mixed_hits")
  expect_equal(result$ids, "20")
  expect_equal(check(blast_test_hits(c("101", "999")))$diagnostics$status, "taxonomy_unavailable")
  expect_equal(check(blast_test_hits(c("1", "1")), list(min_taxa = 1L))$diagnostics$status, "insufficient_resolution")
  one_taxon <- blast_test_hits(c("101", "101"))
  expect_equal(check(one_taxon)$diagnostics$status, "insufficient_records")
  duplicate <- blast_test_hits(); duplicate$sseqid <- "same_record"
  expect_equal(check(duplicate)$diagnostics$status, "insufficient_records")
  weak <- blast_test_hits(); weak$qcovhsp <- 40
  expect_equal(check(weak)$diagnostics$status, "insufficient_quality")
  expect_equal(check(biohelper:::.biohelper_empty_blast_hits())$diagnostics$status, "no_hits")
})

test_that("audit mode flags an outside lineage and unassigns only when requested", {
  paths <- blast_test_paths()
  testthat::local_mocked_bindings(.biohelper_ncbi_id_paths = function(ids, sqlFile) paths[ids])
  matches <- data.frame(feature_id = "asv1", id = "20")
  result <- biohelper:::.biohelper_check_blast(matches, list("20"), "unused", "all", NULL, blast_test_hits(), list())
  expect_equal(result$diagnostics$status, "blast_conflict")
  expect_equal(result$ids, "20")
  result <- biohelper:::.biohelper_check_blast(matches, list("20"), "unused", "all", NULL, blast_test_hits(), list(conflict_action = "unassign"))
  expect_true(is.na(result$ids))
})

test_that("real NCBI ancestry includes intermediate nodes and detects unavailable IDs", {
  database <- ncbi_prunella_database()
  paths <- biohelper:::.biohelper_ncbi_id_paths(c("39176", "175132", "99999999"), database)
  expect_true(paths[["39176"]]$complete)
  expect_true(all(c("39176", "2759", "1") %in% paths[["39176"]]$ids))
  expect_false("175132" %in% paths[["39176"]]$ids)
  expect_false(paths[["99999999"]]$complete)
})

test_that("normalisation integrates precomputed evidence without adding taxonomy columns", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", phylum = "Eukaryota", genus = "Prunella")
  hits <- blast_test_hits(c("175132", "175132"))
  warnings <- character()
  result <- withCallingHandlers(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(min_taxa = 1L), return_details = TRUE
  ), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_equal(result$taxonomy$family, "Prunellidae")
  expect_equal(result$blast$diagnostics$status, "blast_resolved")
  expect_false(any(grepl("Using the first matching ID", warnings, fixed = TRUE)))
  expect_equal(result$blast$hit_source, "precomputed")
  expect_false("status" %in% names(result$taxonomy))
  output <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(min_taxa = 1L)
  ))
  expect_s3_class(output, "data.frame")
  expect_equal(attr(output, "blast_check")$diagnostics$status, "blast_resolved")
})

test_that("BLAST is opt-in and rejects missing sequences or invalid controls", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", genus = "Lamiaceae")
  testthat::local_mocked_bindings(.biohelper_check_blast = function(...) stop("Unexpected BLAST call"))
  expect_no_error(taxo_normalisation(input, database, addExtra = FALSE))
  expect_error(biohelper:::.biohelper_blast_options(list(min_coverage = 101)), "Identity/coverage")
  expect_error(biohelper:::.biohelper_blast_options(list(typo = 5)), "Unknown BLAST")
  expect_error(biohelper:::.biohelper_blast_options(list(poll_seconds = 1)), "at least 60")
  expect_error(biohelper:::.biohelper_blast_sequences(c(asv2 = "ACGT"), "asv1"), "missing")
})

test_that("a real uniquely assigned NCBI family is flagged when BLAST supports birds", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", genus = "Lamiaceae")
  hits <- blast_test_hits(c("175132", "175132"))
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(scope = "all", min_taxa = 1L), return_details = TRUE
  ))
  expect_equal(result$blast$diagnostics$status, "blast_conflict")
  expect_equal(result$taxonomy$family, "Lamiaceae")
  result <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check", blast_hits = hits,
    blast_options = list(scope = "all", min_taxa = 1L, conflict_action = "unassign"), return_details = TRUE
  ))
  expect_true(is.na(result$taxonomy$family))
  expect_equal(result$taxonomy$domain, "Unknown")
})

test_that("phyloseq BLAST diagnostics do not alter abundances or input objects", {
  database <- ncbi_prunella_database()
  abundance <- matrix(10, dimnames = list("asv1", "sample1"))
  taxonomy <- matrix(c("Eukaryota", "Prunella"), nrow = 1,
                     dimnames = list("asv1", c("phylum", "genus")))
  input <- phyloseq::phyloseq(phyloseq::otu_table(abundance, taxa_are_rows = TRUE), phyloseq::tax_table(taxonomy))
  original <- input
  output <- suppressWarnings(taxo_normalisation(
    input, database, addExtra = FALSE, ambiguity = "blast_check",
    blast_hits = blast_test_hits(c("175132", "175132")), blast_options = list(min_taxa = 1L)
  ))
  expect_true(methods::validObject(output))
  expect_identical(phyloseq::otu_table(output), phyloseq::otu_table(input))
  expect_identical(input, original)
  expect_equal(attr(output, "blast_check")$diagnostics$status, "blast_resolved")
})

test_that("search caching deduplicates DNA and remaps hits to each feature", {
  calls <- 0L
  testthat::local_mocked_bindings(.biohelper_run_blast = function(query, output, log, control) {
    calls <<- calls + 1L
    dna <- Biostrings::readDNAStringSet(query)
    expect_length(dna, 1L)
    writeLines(paste(names(dna), 4, 100, "record", 0, "1e-20", 20, 101, "Taxon", "taxon", 100, 100, sep = "\t"), output)
  })
  controls <- biohelper:::.biohelper_blast_options(list(cache_dir = tempfile(), remote = TRUE))
  controls$db <- "core_nt"
  sequences <- c(asv1 = "ACGT", asv2 = "ACGT")
  first <- biohelper:::.biohelper_blast_search(sequences, controls)
  second <- biohelper:::.biohelper_blast_search(sequences, controls)
  expect_equal(calls, 1L)
  expect_equal(first$qseqid, c("asv1", "asv2"))
  expect_identical(first, second)
  controls$max_remote_queries <- 1L
  expect_error(biohelper:::.biohelper_blast_search(c(a = "AAAA", b = "CCCC"), controls), "configured limit")
})

blast_json_fixture <- function(query = "asv1") {
  jsonlite::toJSON(list(BlastOutput2 = list(list(report = list(results = list(search = list(
    query_id = "Query_1", query_title = query, query_len = 100,
    hits = list(list(description = list(list(accession = "record1", taxid = 101, sciname = "Taxon1"),
                                         list(accession = "alias1", taxid = 201, sciname = "Taxon2")),
                     hsps = list(list(query_from = 1, query_to = 80, identity = 80, align_len = 80, evalue = 1e-30, bit_score = 100),
                                 list(query_from = 61, query_to = 100, identity = 38, align_len = 40, evalue = 1e-20, bit_score = 50))))
  )))))), auto_unbox = TRUE)
}

test_that("web parsing preserves merged TaxIDs and single-HSP coverage", {
  testthat::skip_if_not_installed("jsonlite")
  hits <- biohelper:::.biohelper_parse_blast_json(blast_json_fixture(), "asv1")
  expect_equal(hits$pident, c(100, 95))
  expect_equal(hits$qcovhsp, c(80, 40))
  expect_equal(hits$qcovs, c(100, 100))
  expect_equal(hits$staxids, c("101;201", "101;201"))
  expect_error(biohelper:::.biohelper_parse_blast_json(blast_json_fixture(), c("asv1", "asv2")), "missing query")
})

test_that("web requests respect polling and resume saved RIDs without resubmission", {
  testthat::skip_if_not_installed("curl")
  testthat::skip_if_not_installed("jsonlite")
  clock <- 0
  submissions <- 0L
  testthat::local_mocked_bindings(
    .biohelper_blast_now = function() clock,
    .biohelper_blast_wait = function(seconds) { expect_gte(seconds, 60); clock <<- clock + seconds },
    .biohelper_blast_http = function(parameters) {
      if (parameters$CMD == "Put") { submissions <<- submissions + 1L; "RID = RID123\nRTOE = 30" }
      else { expect_equal(parameters$RID, "RID123"); blast_json_fixture() }
    }
  )
  query <- tempfile(fileext = ".fasta"); writeLines(c(">asv1", paste(rep("A", 100), collapse = "")), query)
  output <- tempfile(fileext = ".tsv")
  log <- tempfile()
  controls <- biohelper:::.biohelper_blast_options(list(timeout_seconds = 1))
  controls$db <- "core_nt"
  expect_error(biohelper:::.biohelper_run_remote_blast(query, output, log, controls), "still pending")
  expect_true(file.exists(paste0(output, ".request.rds")))
  controls$timeout_seconds <- 900
  biohelper:::.biohelper_run_remote_blast(query, output, log, controls)
  expect_equal(submissions, 1L)
  expect_equal(nrow(read_blast_taxo_csv(output)), 2L)
})
