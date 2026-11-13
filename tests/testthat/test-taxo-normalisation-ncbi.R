# These integration tests use actual taxonomizr queries against a small SQLite
# database built from authentic NCBI EFetch records (2026-11-13). No lookup
# functions are mocked. See fixtures/ncbi-prunella/provenance.json.
normalise_ncbi_with_warnings <- function(input, database, ...) {
  warnings <- character()
  output <- withCallingHandlers(
    taxo_normalisation(input, sqlFile = database, addExtra = FALSE, ...),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(output = output, warnings = warnings)
}

test_that("the real NCBI fixture contains both Prunella homonyms and full lineages", {
  database <- ncbi_prunella_database()
  expect_warning(
    ids <- taxonomizr::getId("Prunella", database),
    "Multiple taxa ids"
  )
  expect_setequal(strsplit(ids, ",", fixed = TRUE)[[1]], c("39176", "175132"))
  lineage <- taxonomizr::getRawTaxonomy(c("39176", "175132"), database)
  expect_equal(unname(lineage[[1]]["subtribe"]), "Prunellinae")
  expect_equal(unname(lineage[[2]]["superfamily"]), "Passeroidea")
  expect_true("Pentapetalae" %in% unname(lineage[[1]]))

  standard <- taxonomizr::getTaxonomy(c("39176", "175132"), database)
  expect_false(any(c("Prunellinae", "Passeroidea", "Pentapetalae") %in% standard))
})

test_that("real NCBI intermediate ranks resolve misranked input names", {
  database <- ncbi_prunella_database()
  input <- data.frame(
    feature_id = c("plant_subtribe", "bird_superfamily", "plant_clade", "family_in_genus"),
    phylum = NA_character_,
    family = c("Prunellinae", "Passeroidea", "Pentapetalae", NA_character_),
    genus = c("Prunella", "Prunella", "Prunella", "Lamiaceae")
  )
  original <- input
  result <- normalise_ncbi_with_warnings(input, database)
  expect_equal(result$output$feature_id, input$feature_id)
  expect_equal(result$output$family, c("Lamiaceae", "Prunellidae", "Lamiaceae", "Lamiaceae"))
  expect_equal(result$output$genus, c("Prunella", "Prunella", "Prunella", NA_character_))
  expect_equal(result$output$phylum, c("Streptophyta", "Chordata", "Streptophyta", "Streptophyta"))
  expect_false(any(grepl("Using the first matching ID", result$warnings, fixed = TRUE)))
  expect_identical(input, original)

  restricted <- normalise_ncbi_with_warnings(input, database, ranks = c("family", "genus"))
  expect_named(restricted$output, c("feature_id", "family", "genus"))
  expect_equal(restricted$output$family, result$output$family)
  expect_equal(restricted$output$genus, result$output$genus)
})

test_that("real NCBI lookup falls back right to left before resolving ambiguity", {
  database <- ncbi_prunella_database()
  input <- data.frame(
    feature_id = "fallback_to_family_column",
    phylum = NA_character_, order = "Prunellinae", family = "Prunella",
    genus = "NotRecognisedGenus", species = "NotRecognised_name"
  )
  result <- normalise_ncbi_with_warnings(input, database)
  expect_equal(result$output$family, "Lamiaceae")
  expect_equal(result$output$genus, "Prunella")
  expect_true(is.na(result$output$species))
  expect_false(any(grepl("Using the first matching ID", result$warnings, fixed = TRUE)))

  # The matched name in the family column must not be reused as its own
  # supporting parent after unsuccessful genus/species lookups.
  unsupported <- input
  unsupported$order <- "NotRecognisedOrder"
  result <- normalise_ncbi_with_warnings(unsupported, database)
  expect_true(is.na(result$output$genus))
  expect_false(any(grepl("Using the first matching ID", result$warnings, fixed = TRUE)))

  # Unknown intermediate names must not prevent reaching the input phylum.
  input <- data.frame(
    feature_id = "phylum_only", phylum = "Chordata", class = "NotRecognisedClass",
    order = "NotRecognisedOrder", family = "NotRecognisedFamily", genus = "Prunella"
  )
  result <- normalise_ncbi_with_warnings(input, database)
  expect_equal(result$output$family, "Prunellidae")
  expect_equal(result$output$phylum, "Chordata")
  expect_equal(result$output$genus, "Prunella")
})

test_that("real NCBI ambiguity uses all broader names regardless of column labels", {
  database <- ncbi_prunella_database()
  input <- data.frame(
    feature_id = "ambiguous_prunella", kingdom = "Eukaryota",
    phylum = "Eukaryota", genus = "Prunella"
  )
  # Both candidates share Eukaryota, so their names alone remain ambiguous.
  ids <- suppressWarnings(taxonomizr::getId("Prunella", database))
  first_id <- strsplit(ids, ",", fixed = TRUE)[[1]][[1]]
  expected <- taxonomizr::getTaxonomy(first_id, database, desiredTaxa = "family")
  result <- normalise_ncbi_with_warnings(input, database)
  expect_equal(result$output$family, unname(expected[1, 1]))
  expect_equal(result$output$family, "Lamiaceae")
  expect_true(any(grepl(
    "Multiple NCBI IDs.*ambiguous_prunella.*first matching ID", result$warnings
  )))

  # A clue before a column labelled phylum must still resolve the homonym.
  input$kingdom <- "Metazoa"
  input$phylum <- "NotRecognisedPhylum"
  result <- normalise_ncbi_with_warnings(input, database)
  expect_equal(result$output$family, "Prunellidae")
  expect_equal(result$output$phylum, "Chordata")

  # If no broader name supports any candidate, retain the unassigned outcome.
  input$kingdom <- "NotRecognisedKingdom"
  result <- normalise_ncbi_with_warnings(input, database)
  expect_true(is.na(result$output$genus))
  expect_true(is.na(result$output$family))
})

test_that("real NCBI normalisation preserves phyloseq features and abundances", {
  database <- ncbi_prunella_database()
  taxonomy <- matrix(
    c(NA_character_, "Prunellinae", "Prunella"), nrow = 1,
    dimnames = list("asv_plant", c("phylum", "family", "genus"))
  )
  abundance <- matrix(10, dimnames = list("asv_plant", "sample1"))
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(abundance, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy)
  )
  original <- ps
  result <- normalise_ncbi_with_warnings(ps, database)
  expect_s4_class(result$output, "phyloseq")
  expect_equal(phyloseq::taxa_names(result$output), "asv_plant")
  expect_identical(phyloseq::otu_table(result$output), phyloseq::otu_table(ps))
  expect_equal(unname(methods::as(phyloseq::tax_table(result$output), "matrix")[, "family"]), "Lamiaceae")
  expect_identical(ps, original)
})

test_that("real NCBI lookup handles wholly missing or unrecognised taxonomy", {
  database <- ncbi_prunella_database()
  input <- data.frame(
    feature_id = c("missing", "unrecognised"), phylum = NA_character_,
    family = c(NA_character_, "NotRecognisedFamily"),
    genus = c(NA_character_, "NotRecognisedGenus")
  )
  result <- normalise_ncbi_with_warnings(input, database)
  expect_equal(result$output$feature_id, input$feature_id)
  expect_equal(result$output$domain, rep("Unknown", 2))
  expect_true(all(is.na(result$output$family)))
  expect_true(all(is.na(result$output$genus)))
  expect_length(result$warnings, 0L)
})

test_that("large ambiguity warnings stay readable and retain every affected ID", {
  database <- ncbi_prunella_database()
  input <- data.frame(
    feature_id = sprintf("asv_%032d", seq_len(100)),
    phylum = "Eukaryota", genus = "Prunella"
  )
  ambiguity <- NULL
  output <- withCallingHandlers(
    taxo_normalisation(input, sqlFile = database, addExtra = FALSE),
    warning = function(w) {
      if (inherits(w, "biohelper_taxonomy_ambiguity")) ambiguity <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(ambiguity, "biohelper_taxonomy_ambiguity")
  expect_match(conditionMessage(ambiguity), "100 feature\\(s\\)")
  expect_match(conditionMessage(ambiguity), "first 5 shown")
  expect_match(conditionMessage(ambiguity), "Using the first matching ID")
  expect_lt(nchar(conditionMessage(ambiguity)), 1000)
  expect_identical(ambiguity$feature_ids, input$feature_id)
  expect_equal(nrow(output), nrow(input))
  expect_true(all(output$family == "Lamiaceae"))
})
