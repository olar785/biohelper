test_that("representative hits use deterministic E-value/score/identity/coverage/subject ordering", {
  hits <- rbind(
    workflow_hits(subject = "z", evalue = 0, bitscore = 100, identity = 100),
    workflow_hits(subject = "z", evalue = 0, bitscore = 200, identity = 98),
    workflow_hits(subject = "z", evalue = 0, bitscore = 200, identity = 99, coverage = 95),
    workflow_hits(subject = "z", evalue = 0, bitscore = 200, identity = 99, coverage = 100),
    workflow_hits(subject = "a", evalue = 0, bitscore = 200, identity = 99, coverage = 100),
    workflow_hits(subject = "lowest_evalue_loses", evalue = 1e-100, bitscore = 500, identity = 100)
  )
  expect_equal(biohelper:::.biohelper_representative_order(hits), c(5L, 4L, 3L, 2L, 1L, 6L))
  for (indexes in list(seq_len(nrow(hits)), rev(seq_len(nrow(hits))), c(3, 1, 6, 2, 5, 4))) {
    shuffled <- hits[indexes, ]
    expect_equal(shuffled$sseqid[biohelper:::.biohelper_representative_order(shuffled)[[1]]], "a")
  }
  database <- workflow_database()
  result <- lcaPident(workflow_write_hits(hits), tempfile(), taxonomy_db = database)
  expect_equal(result$pident, 99)
  expect_equal(result$qcovs, 100)
  expect_equal(result$evalue, 0)
  expect_equal(result$sseqid, "a")
  expect_equal(result$Percent_Identity, result$pident)
})

test_that("LCA truncates at the broadest conflict and includes every merged-hit TaxID", {
  database <- workflow_database()
  cases <- list(
    same_family = workflow_hits(ids = c("x", "x"), taxids = c("9", "11"), subject = c("a", "b")),
    different_kingdom = workflow_hits(ids = c("x", "x"), taxids = c("9", "15"), subject = c("a", "b")),
    merged_reference = workflow_hits(ids = "x", taxids = "9;11")
  )
  family <- workflow_taxonomy(lcaPident(workflow_write_hits(cases$same_family), tempfile(), taxonomy_db = database))
  expect_equal(family$Family, "Hominidae")
  expect_true(is.na(family$Genus))
  expect_true(is.na(family$Species))
  cross <- workflow_taxonomy(lcaPident(workflow_write_hits(cases$different_kingdom), tempfile(), taxonomy_db = database))
  expect_equal(cross$Domain, "Eukaryota")
  expect_true(all(is.na(cross[, biohelper:::.biohelper_taxonomy_output_ranks()[-1]])))
  merged <- workflow_taxonomy(lcaPident(workflow_write_hits(cases$merged_reference), tempfile(), taxonomy_db = database))
  expect_identical(merged, family)
})

test_that("before, after and no preserve their distinct identity-threshold meanings", {
  database <- workflow_database()
  hits <- workflow_hits(ids = c("x", "x"), taxids = "9", identity = c(90, 98),
                        evalue = c(1e-40, 1e-30), subject = c("low_identity", "high_identity"))
  file <- workflow_write_hits(hits)
  before <- workflow_taxonomy(lcaPident(file, tempfile(), taxonomy_db = database, pident = "before"))
  after <- workflow_taxonomy(lcaPident(file, tempfile(), taxonomy_db = database, pident = "after"))
  no <- workflow_taxonomy(lcaPident(file, tempfile(), taxonomy_db = database, pident = "no", minSim = 100))
  expect_equal(before$Species, "Homo sapiens")
  expect_equal(after$Family, "Hominidae")
  expect_true(is.na(after$Genus))
  expect_true(is.na(after$Species))
  expect_equal(no$Species, "Homo sapiens")
  expect_identical(before, no)
})

test_that("default rank thresholds and coverage boundaries are unchanged", {
  database <- workflow_database()
  thresholds <- c(Kingdom = 71, Phylum = 79, Class = 81, Order = 83, Family = 87, Genus = 95, Species = 97)
  identity <- sort(unique(c(thresholds, thresholds - 0.1, 100)))
  hits <- workflow_hits(ids = paste0("x", seq_along(identity)), identity = identity, coverage = 80)
  file <- workflow_write_hits(hits)
  for (mode in c("before", "after")) {
    result <- lcaPident(file, tempfile(), taxonomy_db = database, pident = mode)
    taxonomy <- workflow_taxonomy(result)
    original_identity <- identity[match(taxonomy$ASV, hits$qseqid)]
    for (rank in names(thresholds)) expect_equal(!is.na(taxonomy[[rank]]), original_identity >= thresholds[[rank]])
  }
  result <- lcaPident(file, tempfile(), taxonomy_db = database, minCov = 80.1)
  expect_equal(result$ASVs, sort(hits$qseqid))
  expect_true(all(is.na(result$pident)))
  expect_true(all(is.na(workflow_taxonomy(result)$Species)))
})

test_that("unresolved IDs do not hide known taxonomy or infrastructure failures", {
  database <- workflow_database()
  hits <- workflow_hits(ids = c("x", "x"), taxids = c("99999999", "9"),
                        evalue = c(0, 1e-30), subject = c("unknown", "known"))
  result <- lcaPident(workflow_write_hits(hits), tempfile(), taxonomy_db = database)
  expect_equal(workflow_taxonomy(result)$Species, "Homo sapiens")
  expect_equal(result$sseqid, "unknown")
  expect_equal(attr(result, "provenance")$unresolved_taxids, 1L)
  expect_equal(attr(result, "provenance")$resolved_taxids, 1L)
  missing <- lcaPident(workflow_write_hits(workflow_hits(taxids = "N/A")), tempfile(), taxonomy_db = database)
  expect_true(is.na(workflow_taxonomy(missing)$Species))
  expect_equal(missing$pident, 99)

  connection <- DBI::dbConnect(RSQLite::SQLite(), database)
  DBI::dbExecute(connection, "UPDATE nodes SET parent=99999999 WHERE id=8")
  DBI::dbDisconnect(connection)
  expect_error(lcaPident(workflow_write_hits(workflow_hits()), tempfile(), taxonomy_db = database), "Corrupt NCBI taxonomy: missing ancestor")
})

test_that("local lookups support taxizedb and taxonomizr schemas and authentic NCBI data", {
  compact <- workflow_database()
  legacy <- workflow_database("taxizedb")
  a <- biohelper:::.biohelper_lca_lookup(c("9", "11", "999"), biohelper:::.biohelper_validate_ncbi_db(compact))
  b <- biohelper:::.biohelper_lca_lookup(c("9", "11", "999"), biohelper:::.biohelper_validate_ncbi_db(legacy))
  expect_identical(a, b)
  authentic <- ncbi_prunella_database()
  result <- lcaPident(workflow_write_hits(workflow_hits(ids = c("x", "x"), taxids = c("39176", "175132"))),
                      tempfile(), taxonomy_db = authentic)
  expect_equal(workflow_taxonomy(result)$Domain, "Eukaryota")
  expect_true(is.na(workflow_taxonomy(result)$Kingdom))
})

test_that("empty and fully filtered BLAST files retain a consistent output schema", {
  database <- workflow_database()
  empty <- lcaPident(workflow_write_hits(workflow_hits()[FALSE, ]), tempfile(), taxonomy_db = database)
  expect_equal(nrow(empty), 0L)
  expect_type(empty$pident, "double")
  expect_type(empty$evalue, "double")
  hits <- workflow_hits(ids = c("002", "001"))
  hits$sscinames <- "uncultured organism"
  filtered <- lcaPident(workflow_write_hits(hits), tempfile(), taxonomy_db = database)
  expect_identical(filtered$ASVs, c("001", "002"))
  expect_true(all(filtered$taxonomy == "Unknown"))
  expect_true(all(is.na(filtered$qcovs)))
})
