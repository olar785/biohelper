merge_ambiguous_input <- function(ids = "asv1") {
  data.frame(feature_id = ids, broad_clue = "Eukaryota", finest_clue = "Prunella")
}

test_that("merging uses sequence evidence once per input and shares cached local searches", {
  database <- ncbi_prunella_database()
  input <- merge_ambiguous_input(c("bird", "plant"))
  original <- input
  other <- input[2:1, ]
  names(other) <- c("ASV", "misleading_rank", "another_rank")
  metadata <- input
  metadata$pident <- 99
  metadata$qcovs <- 100
  metadata$evalue <- 1e-30
  metadata$assignment_method <- "megablast"
  metadata$support_method <- "blastn"
  searches <- 0L
  observed <- NULL
  testthat::local_mocked_bindings(.biohelper_run_blast = function(query, output, log, control) {
    searches <<- searches + 1L
    observed <<- control
    dna <- Biostrings::readDNAStringSet(query)
    records <- lapply(seq_along(dna), function(i) {
      id <- if (as.character(dna[[i]]) == "ACGTACGT") "175132" else "39176"
      data.frame(qseqid = names(dna)[[i]], qlen = 8, pident = 99,
        sseqid = c("reference_a", "reference_b"), sgi = 0, evalue = 1e-30,
        bitscore = 100, staxids = id, sscinames = "Prunella", sblastnames = "test",
        qcovs = 100, qcovhsp = 100)
    })
    utils::write.table(do.call(rbind, records), output, quote = FALSE,
      row.names = FALSE, col.names = FALSE, sep = "\t")
  })
  cache <- tempfile("merge-cache-")
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  result <- suppressWarnings(taxo_merge(
    list(input, other, metadata), database, addExtra = FALSE,
    ambiguity = "blast_check", sequences = c(bird = "ACGTACGT", plant = "TTTTTTTT"),
    db = "/database with spaces/nt", blastapp_path = "/tools/blastn",
    blast_options = list(fallback = "lca", min_taxa = 1L, nthreads = 3L, cache_dir = cache)
  ))
  expect_identical(input, original)
  expect_equal(searches, 1L)
  expect_false(observed$remote)
  expect_equal(observed$db, "/database with spaces/nt")
  expect_equal(observed$blastapp_path, "/tools/blastn")
  expect_equal(observed$nthreads, 3L)
  expect_equal(result$family[match(c("bird", "plant"), result$feature_id)], c("Prunellidae", "Lamiaceae"))
  expect_named(result, c("feature_id", tolower(biohelper:::.biohelper_taxonomy_output_ranks())))
  reports <- attr(result, "normalisation")
  expect_length(reports, 3L)
  expect_true(all(vapply(reports, function(x) all(x$blast$diagnostics$status == "blast_resolved"), logical(1))))
})

test_that("missing sequences and insufficient hits retain LCA through the merge", {
  database <- ncbi_prunella_database()
  input <- merge_ambiguous_input(c("present", "absent"))
  searches <- 0L
  testthat::local_mocked_bindings(.biohelper_blast_search = function(sequences, control) {
    searches <<- searches + 1L
    expect_equal(names(sequences), "present")
    # One accession cannot pass the default evidence screening.
    data.frame(qseqid = "present", sseqid = "one_reference", staxids = "175132",
      pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
  })
  fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(fasta), add = TRUE)
  writeLines(c(">present description", "ACGT"), fasta)
  result <- suppressWarnings(taxo_merge(
    list(silva = input, pr2 = input[2:1, ]), database, addExtra = FALSE,
    ambiguity = "blast_check", sequences = fasta, db = "/databases/nt",
    blast_options = list(fallback = "lca")
  ))
  expect_equal(searches, 1L)
  expect_setequal(result$feature_id, c("present", "absent"))
  expect_equal(result$domain, c("Eukaryota", "Eukaryota"))
  expect_true(all(is.na(result$family)))
  expect_true(all(is.na(result$genus)))
  reports <- attr(result, "normalisation")
  expect_named(reports, c("silva", "pr2"))
  expect_equal(reports$silva$blast$diagnostics$status, c("insufficient_records", "missing_sequence"))
  expect_true(all(reports$silva$ambiguity$action == "lca"))
})

test_that("default merging stays opt-in and non-BLAST ambiguity strategies are forwarded", {
  database <- ncbi_prunella_database()
  input <- merge_ambiguous_input()
  testthat::local_mocked_bindings(.biohelper_check_blast = function(...) stop("Unexpected BLAST call"))
  default <- suppressWarnings(taxo_merge(list(input), database, addExtra = FALSE))
  expect_equal(default$family, "Lamiaceae")
  expect_null(attr(default, "normalisation"))
  lca <- suppressWarnings(taxo_merge(list(input), database, addExtra = FALSE, ambiguity = "lca"))
  expect_equal(lca$domain, "Eukaryota")
  expect_true(is.na(lca$genus))
  unassigned <- suppressWarnings(taxo_merge(
    list(input), database, addExtra = FALSE, ambiguity = "unassigned"
  ))
  expect_equal(nrow(unassigned), 0L)
  expect_named(unassigned, names(lca))
  expect_equal(attr(unassigned, "normalisation")[[1]]$ambiguity$action, "unassigned")
})

test_that("merge validates options early and propagates BLAST execution failures", {
  testthat::local_mocked_bindings(taxo_normalisation = function(...) stop("Unexpected normalisation"))
  expect_error(taxo_merge(list(), "unused"), "non-empty list")
  expect_error(taxo_merge(list(merge_ambiguous_input()), "unused", ranks = "typo"), "recognised")
  expect_error(taxo_merge(list(merge_ambiguous_input()), "unused", ambiguity = "typo"), "arg")
  expect_error(taxo_merge(list(merge_ambiguous_input()), "unused", ambiguity = "blast_check",
    blast_options = list(scope = "typo")), "Invalid BLAST")
})

test_that("a failed BLAST search stops merging rather than silently using LCA", {
  database <- ncbi_prunella_database()
  testthat::local_mocked_bindings(.biohelper_blast_search = function(...) stop("Local search failed"))
  expect_error(suppressWarnings(taxo_merge(
    list(merge_ambiguous_input()), database, addExtra = FALSE,
    ambiguity = "blast_check", sequences = c(asv1 = "ACGT"), db = "/databases/nt",
    blast_options = list(fallback = "lca")
  )), "Local search failed")
})

test_that("merge forwards the requested species formatting", {
  observed <- NULL
  testthat::local_mocked_bindings(taxo_normalisation = function(obj, spnc, ...) {
    observed <<- spnc
    data.frame(feature_id = obj[[1]], domain = "Eukaryota", genus = "Prunella")
  })
  taxo_merge(list(merge_ambiguous_input()), "unused", spnc = TRUE)
  expect_true(observed)
})
