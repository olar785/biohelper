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

worms_like_table <- function() {
  data.frame(
    scientificname = "Salmo salar",
    rank = "Species",
    AphiaID = "127186",
    valid_AphiaID = "127186",
    valid_name = "Salmo salar",
    isMarine = TRUE,
    isFreshwater = TRUE,
    stringsAsFactors = FALSE
  )
}

obis_like_table <- function() {
  data.frame(
    scientificName = "Salmo salar",
    taxonRank = "Species",
    AphiaID = "127186",
    occurrenceID = "obis-1",
    decimalLatitude = -36.5,
    decimalLongitude = 174.8,
    locality = "Hauraki Gulf",
    country = "New Zealand",
    basisOfRecord = "HumanObservation",
    individualCount = 3,
    stringsAsFactors = FALSE
  )
}

gbif_like_table <- function() {
  data.frame(
    scientificName = "Salmo salar",
    taxonRank = "Species",
    taxonKey = "5204019",
    gbifID = "gbif-1",
    decimalLatitude = -36.5,
    decimalLongitude = 174.8,
    locality = "Hauraki Gulf",
    country = "New Zealand",
    basisOfRecord = "PRESERVED_SPECIMEN",
    individualCount = 1,
    stringsAsFactors = FALSE
  )
}

bold_like_table <- function() {
  data.frame(
    species_name = "Salmo salar",
    genus_name = "Salmo",
    family_name = "Salmonidae",
    processid = "BOLD-1",
    country = "New Zealand",
    province_state = "Auckland",
    exactsite = "Hauraki Gulf",
    lat = -36.5,
    lon = 174.8,
    stringsAsFactors = FALSE
  )
}

test_that("empty build_taxon_evidence returns an empty standard evidence table", {
  evidence <- build_taxon_evidence()

  expect_s3_class(evidence, "data.frame")
  expect_equal(nrow(evidence), 0)
  expect_identical(colnames(evidence), biohelper:::.taxon_evidence_standard_columns())
})

test_that("user_evidence is accepted and returned", {
  user_evidence <- build_taxon_evidence_user_evidence()

  evidence <- build_taxon_evidence(user_evidence = user_evidence)

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$source, "literature")
  expect_equal(evidence$doi, "10.0000/example")
  expect_true(all(biohelper:::.taxon_evidence_standard_columns() %in% colnames(evidence)))
})

test_that("invalid user_evidence errors clearly", {
  user_evidence <- build_taxon_evidence_user_evidence()
  user_evidence$reference <- NULL

  expect_error(
    build_taxon_evidence(user_evidence = user_evidence),
    "taxon_evidence.*missing required columns: reference"
  )
})

test_that("sources must be a named list or error clearly", {
  expect_error(
    build_taxon_evidence(sources = list(worms_like_table())),
    "`sources` must be a named list"
  )

  expect_error(
    build_taxon_evidence(sources = list(worms = "not a data.frame")),
    "`sources\\$worms` must be a data.frame"
  )
})

test_that("worms-like table is standardised correctly", {
  evidence <- build_taxon_evidence(sources = list(worms = worms_like_table()))

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$taxon_rank, "Species")
  expect_equal(evidence$source, "worms")
  expect_equal(evidence$evidence_type, "taxonomic_database")
  expect_equal(evidence$accepted_name, "Salmo salar")
  expect_equal(evidence$source_taxon_id, "127186")
  expect_equal(evidence$environment, "marine; freshwater")
  expect_equal(evidence$reference, "WoRMS")
  expect_match(evidence$reference_url, "127186")
})

test_that("source_type can describe a single source data.frame", {
  evidence <- build_taxon_evidence(
    sources = worms_like_table(),
    source_type = "worms"
  )

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$source, "worms")
  expect_equal(evidence$evidence_type, "taxonomic_database")
})

test_that("obis-like table is standardised correctly", {
  evidence <- build_taxon_evidence(sources = list(obis = obis_like_table()))

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$source, "obis")
  expect_equal(evidence$evidence_type, "occurrence")
  expect_equal(evidence$source_taxon_id, "127186")
  expect_equal(evidence$source_record_id, "obis-1")
  expect_equal(evidence$environment, "marine")
  expect_equal(evidence$region, "New Zealand")
  expect_equal(evidence$locality, "Hauraki Gulf")
  expect_equal(evidence$occurrence_count, "3")
  expect_equal(evidence$reference, "OBIS")
})

test_that("gbif-like table is standardised correctly", {
  evidence <- build_taxon_evidence(sources = list(gbif = gbif_like_table()))

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$source, "gbif")
  expect_equal(evidence$evidence_type, "occurrence")
  expect_equal(evidence$source_taxon_id, "5204019")
  expect_equal(evidence$source_record_id, "gbif-1")
  expect_equal(evidence$region, "New Zealand")
  expect_equal(evidence$locality, "Hauraki Gulf")
  expect_equal(evidence$basis_of_record, "PRESERVED_SPECIMEN")
  expect_equal(evidence$occurrence_count, "1")
  expect_equal(evidence$reference, "GBIF")
})

test_that("bold-like table is standardised correctly", {
  evidence <- build_taxon_evidence(sources = list(bold = bold_like_table()))

  expect_equal(nrow(evidence), 1)
  expect_equal(evidence$taxon_name, "Salmo salar")
  expect_equal(evidence$taxon_rank, "species")
  expect_equal(evidence$source, "bold")
  expect_equal(evidence$evidence_type, "barcode_or_specimen")
  expect_equal(evidence$source_record_id, "BOLD-1")
  expect_equal(evidence$region, "New Zealand; Auckland")
  expect_equal(evidence$locality, "Hauraki Gulf")
  expect_equal(evidence$decimal_latitude, "-36.5")
  expect_equal(evidence$decimal_longitude, "174.8")
  expect_equal(evidence$reference, "BOLD")
})

test_that("combined sources return one evidence table", {
  evidence <- build_taxon_evidence(
    sources = list(
      worms = worms_like_table(),
      obis = obis_like_table(),
      gbif = gbif_like_table(),
      bold = bold_like_table()
    ),
    user_evidence = build_taxon_evidence_user_evidence()
  )

  expect_equal(nrow(evidence), 5)
  expect_equal(
    evidence$source,
    c("worms", "obis", "gbif", "bold", "literature")
  )
  expect_true(all(biohelper:::.taxon_evidence_standard_columns() %in% colnames(evidence)))
})

test_that("no external API calls are made", {
  evidence <- build_taxon_evidence(
    sources = list(worms = worms_like_table())
  )

  expect_equal(evidence$source, "worms")
  expect_equal(evidence$evidence_type, "taxonomic_database")
  expect_false("query_live" %in% names(formals(build_taxon_evidence)))
})
