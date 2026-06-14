fake_worms_db_taxonomy <- function() {
  data.frame(
    AphiaID = c("127186", "126436", "127186"),
    scientificName = c("Salmo salar", "Gadus morhua", "Salmo salar"),
    taxonRank = c("Species", "Species", "Species"),
    stringsAsFactors = FALSE
  )
}

fake_worms_db_evidence <- function(taxa, by = "aphia_id", checked_at = "2026-01-01") {
  evidence <- biohelper:::.empty_taxon_evidence_table(length(taxa))
  evidence$taxon_name <- if (identical(by, "name")) taxa else paste0("taxon_", taxa)
  evidence$taxon_rank <- "Species"
  evidence$source <- "worms"
  evidence$evidence_type <- "taxonomic_environment_database"
  evidence$evidence_summary <- paste("Mock WoRMS evidence for", taxa)
  evidence$reference <- "WoRMS"
  evidence$accepted_name <- evidence$taxon_name
  evidence$accepted_rank <- "Species"
  evidence$source_taxon_id <- if (identical(by, "aphia_id")) taxa else paste0("id_", seq_along(taxa))
  evidence$source_record_id <- evidence$source_taxon_id
  evidence$environment <- "marine"
  evidence$reference_url <- paste0(
    "https://www.marinespecies.org/aphia.php?p=taxdetails&id=",
    evidence$source_taxon_id
  )
  evidence$checked_at <- as.character(checked_at)
  evidence
}

test_that("build_worms_evidence_db reads data.frame input", {
  output_path <- tempfile(fileext = ".rds")
  calls <- list()

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      calls[[length(calls) + 1]] <<- list(taxa = taxa, by = by)
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    batch_size = 100,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(nrow(evidence), 2)
  expect_equal(calls[[1]]$by, "aphia_id")
  expect_equal(calls[[1]]$taxa, c("127186", "126436"))
})

test_that("build_worms_evidence_db reads csv and rds input", {
  csv_path <- tempfile(fileext = ".csv")
  rds_path <- tempfile(fileext = ".rds")
  csv_output <- tempfile(fileext = ".rds")
  rds_output <- tempfile(fileext = ".rds")
  utils::write.csv(fake_worms_db_taxonomy(), csv_path, row.names = FALSE)
  saveRDS(fake_worms_db_taxonomy(), rds_path)

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  csv_evidence <- build_worms_evidence_db(
    csv_path,
    output_path = csv_output,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )
  rds_evidence <- build_worms_evidence_db(
    rds_path,
    output_path = rds_output,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(nrow(csv_evidence), 2)
  expect_equal(nrow(rds_evidence), 2)
  expect_true(file.exists(csv_output))
  expect_true(file.exists(rds_output))
})

test_that("build_worms_evidence_db detects AphiaID column", {
  output_path <- tempfile(fileext = ".rds")
  query_mode <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      query_mode <<- by
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(query_mode, "aphia_id")
})

test_that("build_worms_evidence_db falls back to name column", {
  output_path <- tempfile(fileext = ".rds")
  taxonomy <- data.frame(
    scientificName = c("Salmo salar", "Gadus morhua"),
    taxonRank = c("Species", "Species"),
    stringsAsFactors = FALSE
  )
  query_mode <- NULL
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      query_mode <<- by
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  build_worms_evidence_db(
    taxonomy,
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(query_mode, "name")
  expect_equal(queried_taxa, c("Salmo salar", "Gadus morhua"))
})

test_that("build_worms_evidence_db deduplicates taxa", {
  output_path <- tempfile(fileext = ".rds")
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(queried_taxa, c("127186", "126436"))
})

test_that("build_worms_evidence_db writes output_path", {
  output_path <- tempfile(fileext = ".rds")

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_true(file.exists(output_path))
  expect_equal(readRDS(output_path), evidence)
})

test_that("resume = TRUE skips already fetched taxa", {
  output_path <- tempfile(fileext = ".rds")
  existing <- fake_worms_db_evidence("127186", by = "aphia_id")
  saveRDS(existing, output_path)
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(queried_taxa, "126436")
  expect_equal(nrow(evidence), 2)
  expect_equal(evidence$source_taxon_id, c("127186", "126436"))
})

test_that("overwrite = FALSE errors when output exists and resume = FALSE", {
  output_path <- tempfile(fileext = ".rds")
  saveRDS(fake_worms_db_evidence("127186", by = "aphia_id"), output_path)

  expect_error(
    build_worms_evidence_db(
      fake_worms_db_taxonomy(),
      output_path = output_path,
      resume = FALSE,
      overwrite = FALSE,
      sleep = 0,
      verbose = FALSE
    ),
    "`output_path` already exists"
  )
})

test_that("overwrite = TRUE rebuilds", {
  output_path <- tempfile(fileext = ".rds")
  saveRDS(fake_worms_db_evidence("old", by = "aphia_id"), output_path)
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    overwrite = TRUE,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(queried_taxa, c("127186", "126436"))
  expect_equal(nrow(evidence), 2)
  expect_false("old" %in% evidence$source_taxon_id)
})

test_that("build_worms_evidence_db output passes validate_taxon_evidence", {
  output_path <- tempfile(fileext = ".rds")

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_no_error(biohelper:::validate_taxon_evidence(evidence))
})
