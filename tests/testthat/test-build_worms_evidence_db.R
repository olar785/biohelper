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
  evidence$source_taxon_id <- if (identical(by, "aphia_id")) {
    as.character(taxa)
  } else {
    paste0("id_", seq_along(taxa))
  }
  evidence$source_record_id <- evidence$source_taxon_id
  evidence$environment <- "marine"
  evidence$reference_url <- paste0(
    "https://www.marinespecies.org/aphia.php?p=taxdetails&id=",
    evidence$source_taxon_id
  )
  evidence$checked_at <- as.character(checked_at)
  evidence$fetch_extra_column <- paste0("extra_", as.character(taxa))
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
  expect_equal(calls[[1]]$taxa, c(127186, 126436))
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

test_that("build_worms_evidence_db uses AphiaID when both AphiaID and name columns exist", {
  output_path <- tempfile(fileext = ".rds")
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
    fake_worms_db_taxonomy(),
    output_path = output_path,
    name_col = "scientificName",
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(query_mode, "aphia_id")
  expect_equal(queried_taxa, c(127186, 126436))
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

  expect_equal(queried_taxa, c(127186, 126436))
})

test_that("build_worms_evidence_db accepts CSV-derived character AphiaIDs", {
  csv_path <- tempfile(fileext = ".csv")
  output_path <- tempfile(fileext = ".rds")
  writeLines(
    c(
      "AphiaID,scientificName",
      "1054700,Taxon one",
      "not_an_id,Taxon two",
      " 127186,Taxon three"
    ),
    csv_path
  )
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  expect_warning(
    evidence <- build_worms_evidence_db(
      csv_path,
      output_path = output_path,
      sleep = 0,
      checked_at = "2026-01-01",
      verbose = FALSE
    ),
    "Dropped 1 missing or invalid AphiaID value"
  )

  expect_equal(queried_taxa, c(1054700, 127186))
  expect_equal(evidence$source_taxon_id, c("1054700", "127186"))
})

test_that("invalid AphiaIDs are dropped with a clear warning", {
  output_path <- tempfile(fileext = ".rds")
  taxonomy <- data.frame(
    AphiaID = c("bad", NA, "", "126436"),
    scientificName = c("Bad one", "Missing one", "Empty one", "Gadus morhua"),
    stringsAsFactors = FALSE
  )
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  expect_warning(
    evidence <- build_worms_evidence_db(
      taxonomy,
      output_path = output_path,
      sleep = 0,
      checked_at = "2026-01-01",
      verbose = FALSE
    ),
    "Dropped 3 missing or invalid AphiaID value"
  )

  expect_equal(queried_taxa, 126436)
  expect_equal(evidence$source_taxon_id, "126436")
})

test_that("all invalid AphiaIDs error clearly", {
  output_path <- tempfile(fileext = ".rds")
  taxonomy <- data.frame(
    AphiaID = c("bad", NA, ""),
    scientificName = c("Bad one", "Missing one", "Empty one"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    expect_error(
      build_worms_evidence_db(
        taxonomy,
        output_path = output_path,
        sleep = 0,
        checked_at = "2026-01-01",
        verbose = FALSE
      ),
      "All AphiaID values.*missing or invalid.*Check `id_col`"
    ),
    "Dropped 3 missing or invalid AphiaID value"
  )
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

test_that("build_worms_evidence_db saves after each batch", {
  output_path <- tempfile(fileext = ".rds")
  calls <- 0

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      calls <<- calls + 1
      if (calls == 2) {
        expect_true(file.exists(output_path))
        expect_equal(nrow(readRDS(output_path)), 1)
      }
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    batch_size = 1,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(calls, 2)
  expect_equal(nrow(evidence), 2)
  expect_equal(nrow(readRDS(output_path)), 2)
})

test_that("resume = TRUE skips already fetched AphiaIDs across character and numeric IDs", {
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

  expect_equal(queried_taxa, 126436)
  expect_equal(nrow(evidence), 2)
  expect_equal(evidence$source_taxon_id, c("127186", "126436"))
})

test_that("resume = TRUE skips already fetched names when using name fallback", {
  output_path <- tempfile(fileext = ".rds")
  taxonomy <- data.frame(
    scientificName = c("Salmo salar", "Gadus morhua", "Salmo salar"),
    taxonRank = c("Species", "Species", "Genus"),
    stringsAsFactors = FALSE
  )
  existing <- fake_worms_db_evidence("Salmo salar", by = "name")
  saveRDS(existing, output_path)
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    taxonomy,
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(queried_taxa, "Gadus morhua")
  expect_equal(nrow(evidence), 2)
  expect_equal(evidence$taxon_name, c("Salmo salar", "Gadus morhua"))
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

  expect_equal(queried_taxa, c(127186, 126436))
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

test_that("build_worms_evidence_db does not filter by rank", {
  output_path <- tempfile(fileext = ".rds")
  taxonomy <- data.frame(
    AphiaID = c("1", "2", "3"),
    scientificName = c("Species one", "Genus two", "Family three"),
    taxonRank = c("Species", "Genus", "Family"),
    stringsAsFactors = FALSE
  )
  queried_taxa <- NULL

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      queried_taxa <<- taxa
      fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
    }
  )

  evidence <- build_worms_evidence_db(
    taxonomy,
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_equal(queried_taxa, c(1, 2, 3))
  expect_equal(nrow(evidence), 3)
})

test_that("build_worms_evidence_db preserves full fetch output columns", {
  output_path <- tempfile(fileext = ".rds")

  testthat::local_mocked_bindings(
    fetch_worms_evidence = function(taxa, by, marine_only, cache_path, sleep, checked_at, verbose) {
      out <- fake_worms_db_evidence(taxa, by = by, checked_at = checked_at)
      out$reference <- paste0("Reference for ", taxa)
      out$reference_url <- paste0("https://example.org/aphia/", taxa)
      out
    }
  )

  evidence <- build_worms_evidence_db(
    fake_worms_db_taxonomy(),
    output_path = output_path,
    sleep = 0,
    checked_at = "2026-01-01",
    verbose = FALSE
  )

  expect_true("fetch_extra_column" %in% colnames(evidence))
  expect_equal(evidence$reference, paste0("Reference for ", c(127186, 126436)))
  expect_equal(evidence$reference_url, paste0("https://example.org/aphia/", c(127186, 126436)))
  expect_equal(readRDS(output_path)$fetch_extra_column, evidence$fetch_extra_column)
})
