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

mock_worms_records_for_ids <- function(ids) {
  do.call(rbind, lapply(ids, function(id) {
    mock_worms_record(
      scientificname = paste0("taxon_", id),
      AphiaID = as.character(id),
      valid_AphiaID = as.character(id),
      valid_name = paste0("taxon_", id)
    )
  }))
}

strict_worms_record_mock <- function(expected_id) {
  force(expected_id)
  calls <- 0
  list(
    get_calls = function() calls,
    record = function(...) {
      calls <<- calls + 1
      args <- list(...)
      expect_true("id" %in% names(args))
      expect_false(any(c("aphia_id", "AphiaID", "taxa", "x") %in% names(args)))
      expect_true(is.numeric(args$id) || is.integer(args$id))
      expect_equal(args$id, expected_id)
      mock_worms_records_for_ids(args$id)
    }
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
  strict_mock <- strict_worms_record_mock(127186)

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = strict_mock$record
  )

  evidence <- fetch_worms_evidence(
    taxa = 127186,
    by = "aphia_id",
    verbose = FALSE
  )

  expect_identical(colnames(evidence), biohelper:::.taxon_evidence_standard_columns())
  expect_equal(evidence$taxon_name, "taxon_127186")
  expect_equal(evidence$source_taxon_id, "127186")
  expect_equal(evidence$source_record_id, "127186")
  expect_equal(evidence$accepted_name, "taxon_127186")
  expect_equal(evidence$source, "worms")
  expect_s3_class(biohelper:::validate_taxon_evidence(evidence), "data.frame")
  expect_equal(strict_mock$get_calls(), 1)
})

test_that("AphiaID query accepts a single character ID and queries WoRMS with numeric id", {
  strict_mock <- strict_worms_record_mock(1054700)

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = strict_mock$record
  )

  evidence <- fetch_worms_evidence(
    taxa = "1054700",
    by = "aphia_id",
    verbose = FALSE
  )

  expect_equal(evidence$source_taxon_id, "1054700")
  expect_match(evidence$evidence_summary, "Queried taxon '1054700'")
  expect_equal(strict_mock$get_calls(), 1)
})

test_that("AphiaID query accepts a single numeric ID", {
  strict_mock <- strict_worms_record_mock(370514)

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = strict_mock$record
  )

  evidence <- fetch_worms_evidence(
    taxa = 370514,
    by = "aphia_id",
    verbose = FALSE
  )

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$source_taxon_id, "370514")
  expect_equal(strict_mock$get_calls(), 1)
})

test_that("AphiaID query accepts a numeric vector and queries WoRMS once", {
  strict_mock <- strict_worms_record_mock(c(370514, 1054700))

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = strict_mock$record
  )

  evidence <- fetch_worms_evidence(
    taxa = c(370514, 1054700),
    by = "aphia_id",
    verbose = FALSE
  )

  expect_equal(evidence$source_taxon_id, c("370514", "1054700"))
  expect_equal(evidence$taxon_name, c("taxon_370514", "taxon_1054700"))
  expect_equal(strict_mock$get_calls(), 1)
})

test_that("AphiaID query accepts a character vector and queries WoRMS once with numeric id", {
  strict_mock <- strict_worms_record_mock(c(370514, 1054700))

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = strict_mock$record
  )

  evidence <- fetch_worms_evidence(
    taxa = c("370514", "1054700"),
    by = "aphia_id",
    verbose = FALSE
  )

  expect_equal(evidence$source_taxon_id, c("370514", "1054700"))
  expect_equal(evidence$taxon_name, c("taxon_370514", "taxon_1054700"))
  expect_equal(strict_mock$get_calls(), 1)
})

test_that("AphiaID vector query returns no-match rows for missing returned records", {
  calls <- 0

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = function(...) {
      calls <<- calls + 1
      args <- list(...)
      expect_true("id" %in% names(args))
      expect_equal(args$id, c(370514, 1054700))
      mock_worms_records_for_ids(370514)
    }
  )

  evidence <- fetch_worms_evidence(
    taxa = c(370514, 1054700),
    by = "aphia_id",
    verbose = FALSE
  )

  expect_equal(calls, 1)
  expect_equal(evidence$source_taxon_id, c("370514", "1054700"))
  expect_equal(evidence$taxon_name, c("taxon_370514", "1054700"))
  expect_equal(
    evidence$evidence_summary[[2]],
    "No WoRMS record found for queried taxon."
  )
})

test_that("AphiaID query drops invalid IDs with a clear warning", {
  strict_mock <- strict_worms_record_mock(370514)

  testthat::local_mocked_bindings(
    .require_worrms = function() TRUE,
    .call_worrms_wm_record = strict_mock$record
  )

  expect_warning(
    evidence <- fetch_worms_evidence(
      taxa = c("370514", "not_an_id"),
      by = "aphia_id",
      verbose = FALSE
    ),
    "Dropped 1 missing or invalid AphiaID value"
  )

  expect_equal(evidence$source_taxon_id, "370514")
  expect_equal(strict_mock$get_calls(), 1)
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
