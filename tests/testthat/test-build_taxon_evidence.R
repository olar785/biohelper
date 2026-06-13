build_taxon_evidence_taxonomy <- function() {
  data.frame(
    feature_id = c("asv1", "asv2", "asv3"),
    kingdom = c("Animalia", "Animalia", "Animalia"),
    phylum = c("Chordata", "Chordata", "Arthropoda"),
    genus = c("Salmo", "Salmo", "Daphnia"),
    species = c("Salmo salar", "Salmo salar", NA_character_),
    stringsAsFactors = FALSE
  )
}

build_taxon_evidence_user_evidence <- function() {
  data.frame(
    taxon_name = "Salmo salar",
    taxon_rank = "species",
    source = "literature",
    evidence_type = "environment",
    evidence_summary = "Anadromous species with marine and freshwater phases.",
    reference = "Example reference",
    environment = "marine; freshwater",
    region = "North Atlantic",
    doi = "10.0000/example",
    stringsAsFactors = FALSE
  )
}

test_that("build_taxon_evidence returns required evidence columns", {
  evidence <- build_taxon_evidence(build_taxon_evidence_taxonomy())

  expect_s3_class(evidence, "data.frame")
  expect_true(all(biohelper:::.taxon_evidence_required_columns() %in% colnames(evidence)))
  expect_true(all(biohelper:::.taxon_evidence_optional_columns() %in% colnames(evidence)))
})

test_that("one placeholder row is created per unique query taxon", {
  evidence <- build_taxon_evidence(build_taxon_evidence_taxonomy())

  expect_equal(nrow(evidence), 2)
  expect_equal(evidence$evidence_type, rep("not_queried", 2))
  expect_equal(evidence$evidence_summary, rep("No evidence fetched yet.", 2))
  expect_true(all(is.na(evidence$reference)))
  expect_false(any(is.na(evidence$checked_at)))
})

test_that("duplicated taxa are deduplicated", {
  evidence <- build_taxon_evidence(build_taxon_evidence_taxonomy())
  taxon_pairs <- paste(evidence$taxon_name, evidence$taxon_rank, sep = "|")

  expect_equal(length(taxon_pairs), length(unique(taxon_pairs)))
  expect_true("Salmo salar|species" %in% taxon_pairs)
  expect_true("Daphnia|genus" %in% taxon_pairs)
})

test_that("user_evidence is accepted and included", {
  user_evidence <- build_taxon_evidence_user_evidence()

  evidence <- build_taxon_evidence(
    build_taxon_evidence_taxonomy(),
    user_evidence = user_evidence
  )

  expect_true(any(evidence$evidence_type == "environment"))
  expect_true(any(evidence$reference == "Example reference", na.rm = TRUE))
  expect_true(any(evidence$doi == "10.0000/example", na.rm = TRUE))
})

test_that("invalid user_evidence errors clearly", {
  user_evidence <- build_taxon_evidence_user_evidence()
  user_evidence$reference <- NULL

  expect_error(
    build_taxon_evidence(
      build_taxon_evidence_taxonomy(),
      user_evidence = user_evidence
    ),
    "taxon_evidence.*missing required columns: reference"
  )
})

test_that("include_empty = FALSE returns only user_evidence if provided", {
  user_evidence <- build_taxon_evidence_user_evidence()

  evidence <- build_taxon_evidence(
    build_taxon_evidence_taxonomy(),
    user_evidence = user_evidence,
    include_empty = FALSE
  )

  expect_equal(nrow(evidence), nrow(user_evidence))
  expect_equal(evidence$taxon_name, user_evidence$taxon_name)
  expect_equal(evidence$evidence_type, user_evidence$evidence_type)
  expect_false(any(evidence$evidence_type == "not_queried"))
})

test_that("include_empty = FALSE and user_evidence = NULL returns empty standard evidence table", {
  evidence <- build_taxon_evidence(
    build_taxon_evidence_taxonomy(),
    include_empty = FALSE
  )

  expect_s3_class(evidence, "data.frame")
  expect_equal(nrow(evidence), 0)
  expect_identical(colnames(evidence), biohelper:::.taxon_evidence_standard_columns())
})

test_that("no external API calls are made", {
  evidence <- build_taxon_evidence(
    build_taxon_evidence_taxonomy(),
    evidence_sources = "local_test_source"
  )

  expect_equal(unique(evidence$source), "local_test_source")
  expect_true(all(evidence$evidence_type == "not_queried"))
})
