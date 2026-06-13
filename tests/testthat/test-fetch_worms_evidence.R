mock_worms_record <- function(
  scientificname = "Salmo salar",
  rank = "Species",
  AphiaID = "127186",
  valid_AphiaID = "127186",
  valid_name = scientificname,
  status = "accepted",
  isMarine = TRUE,
  isBrackish = FALSE,
  isFreshwater = FALSE,
  isTerrestrial = FALSE,
  citation = "WoRMS Editorial Board",
  url = NA_character_
) {
  data.frame(
    scientificname = scientificname,
    rank = rank,
    AphiaID = AphiaID,
    valid_AphiaID = valid_AphiaID,
    valid_name = valid_name,
    status = status,
    isMarine = isMarine,
    isBrackish = isBrackish,
    isFreshwater = isFreshwater,
    isTerrestrial = isTerrestrial,
    citation = citation,
    url = url,
    stringsAsFactors = FALSE
  )
}

test_that("successful name query returns standard evidence columns", {
  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .worms_records_names = function(taxa, marine_only = FALSE) {
      expect_equal(taxa, "Salmo salar")
      expect_false(marine_only)
      list(mock_worms_record())
    }
  )

  evidence <- fetch_worms_evidence("Salmo salar", verbose = FALSE)

  expect_identical(colnames(evidence), biohelper:::.taxon_evidence_standard_columns())
  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$source, "worms")
  expect_equal(evidence$evidence_type, "taxonomic_environment_database")
  expect_equal(evidence$reference, "WoRMS Editorial Board")
  expect_equal(evidence$source_taxon_id, "127186")
  expect_match(evidence$reference_url, "127186")
  expect_s3_class(biohelper:::validate_taxon_evidence(evidence), "data.frame")
})

test_that("successful AphiaID query returns standard evidence columns", {
  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .worms_record = function(aphia_id) {
      expect_equal(aphia_id, "127186")
      mock_worms_record()
    }
  )

  evidence <- fetch_worms_evidence(
    taxa = 127186,
    by = "aphia_id",
    verbose = FALSE
  )

  expect_identical(colnames(evidence), biohelper:::.taxon_evidence_standard_columns())
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$source_taxon_id, "127186")
  expect_equal(evidence$source_record_id, "127186")
  expect_equal(evidence$accepted_name, "Salmo salar")
  expect_equal(evidence$source, "worms")
  expect_s3_class(biohelper:::validate_taxon_evidence(evidence), "data.frame")
})

test_that("environment is derived from WoRMS environment flags", {
  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .worms_records_names = function(taxa, marine_only = FALSE) {
      list(mock_worms_record(
        isMarine = 1,
        isBrackish = TRUE,
        isFreshwater = 0,
        isTerrestrial = "yes"
      ))
    }
  )

  evidence <- fetch_worms_evidence("Salmo salar", verbose = FALSE)

  expect_equal(evidence$environment, "marine; brackish; terrestrial")
  expect_match(evidence$evidence_summary, "environment flags: marine; brackish; terrestrial")
})

test_that("no match returns a standard evidence row with environment NA", {
  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .worms_records_names = function(taxa, marine_only = FALSE) {
      list(data.frame())
    }
  )

  evidence <- fetch_worms_evidence("Imaginary taxon", verbose = FALSE)

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Imaginary taxon")
  expect_equal(evidence$taxon_rank, "unknown")
  expect_equal(evidence$source, "worms")
  expect_equal(evidence$evidence_type, "taxonomic_environment_database")
  expect_equal(evidence$evidence_summary, "No WoRMS record found for queried taxon.")
  expect_true(is.na(evidence$environment))
  expect_s3_class(biohelper:::validate_taxon_evidence(evidence), "data.frame")
})

test_that("cache_path avoids repeated queries", {
  cache_path <- tempfile(fileext = ".rds")
  query_count <- 0

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .worms_records_names = function(taxa, marine_only = FALSE) {
      query_count <<- query_count + 1
      list(mock_worms_record())
    }
  )

  first <- fetch_worms_evidence(
    "Salmo salar",
    cache_path = cache_path,
    verbose = FALSE
  )
  expect_equal(query_count, 1)

  testthat::local_mocked_bindings(
    .require_worrms = function() {
      stop("worrms should not be required for cached taxa", call. = FALSE)
    },
    .worms_records_names = function(...) {
      stop("WoRMS should not be queried for cached taxa", call. = FALSE)
    }
  )

  second <- fetch_worms_evidence(
    "Salmo salar",
    cache_path = cache_path,
    verbose = FALSE
  )

  expect_equal(first, second)
  expect_equal(query_count, 1)
})

test_that("missing worrms package errors clearly when live query is attempted", {
  testthat::local_mocked_bindings(
    .require_worrms = function() {
      stop(
        "Live WoRMS queries require the `worrms` package. Install it with `install.packages(\"worrms\")`.",
        call. = FALSE
      )
    }
  )

  expect_error(
    fetch_worms_evidence("Salmo salar", verbose = FALSE),
    "worrms.*install.packages"
  )
})

test_that("output passes validate_taxon_evidence", {
  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .worms_records_names = function(taxa, marine_only = FALSE) {
      list(mock_worms_record())
    }
  )

  evidence <- fetch_worms_evidence("Salmo salar", verbose = FALSE)

  expect_no_error(biohelper:::validate_taxon_evidence(evidence))
})
