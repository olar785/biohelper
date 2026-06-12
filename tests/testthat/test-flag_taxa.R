flag_taxa_test_taxonomy <- function() {
  data.frame(
    feature_id = c("asv1", "asv2"),
    kingdom = c("Animalia", "Animalia"),
    phylum = c("Chordata", "Arthropoda"),
    genus = c("Salmo", "Daphnia"),
    species = c("Salmo salar", NA_character_),
    stringsAsFactors = FALSE
  )
}

valid_flag_taxa_result <- function() {
  cbind(
    flag_taxa_test_taxonomy(),
    data.frame(
      query_rank = c("species", "genus"),
      query_name = c("Salmo salar", "Daphnia"),
      lineage = c(
        "kingdom: Animalia; phylum: Chordata; genus: Salmo; species: Salmo salar",
        "kingdom: Animalia; phylum: Arthropoda; genus: Daphnia"
      ),
      expected_environment = c("marine", "marine"),
      expected_environment_status = c("compatible", "mixed_within_rank"),
      expected_habitat = c("estuary", "estuary"),
      expected_habitat_status = c("transient_or_allochthonous_possible", "unknown"),
      expected_region = c("North Atlantic", "North Atlantic"),
      expected_region_status = c("known_in_region", "not_assessed"),
      recommended_action = c("retain", "review"),
      rationale = c("Example rationale 1.", "Example rationale 2."),
      references = c("Example reference 1.", "Example reference 2."),
      stringsAsFactors = FALSE
    )
  )
}

test_that("prompt_only = TRUE returns a character prompt", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    family = "Salmonidae",
    genus = "Salmo",
    species = "Salmo salar"
  )

  prompt <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic"
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
  expect_true(grepl("Salmo salar", prompt, fixed = TRUE))
})

test_that("prompt contains expected context values", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    genus = "Salmo",
    species = "Salmo salar"
  )

  prompt <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "Aotearoa New Zealand"
  )

  expect_true(grepl("marine", prompt, fixed = TRUE))
  expect_true(grepl("estuary", prompt, fixed = TRUE))
  expect_true(grepl("Aotearoa New Zealand", prompt, fixed = TRUE))
})

test_that("prompt contains allowed status values", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    genus = "Salmo",
    species = "Salmo salar"
  )

  prompt <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )

  allowed_values <- c(
    "compatible",
    "incompatible",
    "mixed_within_rank",
    "unknown",
    "insufficient_taxonomic_resolution",
    "transient_or_allochthonous_possible",
    "known_in_region",
    "known_near_region",
    "known_elsewhere_only",
    "restricted_elsewhere",
    "no_distribution_evidence",
    "not_assessed",
    "retain",
    "review",
    "exclude"
  )

  for (allowed_value in allowed_values) {
    expect_true(grepl(allowed_value, prompt, fixed = TRUE))
  }
})

test_that("prompt_only = FALSE errors clearly", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    species = "Salmo salar"
  )

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE
    ),
    "LLM backend.*not implemented"
  )
})

test_that("valid mock_llm_result is returned", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_identical(output, result)
})

test_that("invalid recommended_action errors", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$recommended_action[1] <- "maybe"

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "recommended_action.*invalid values.*maybe"
  )
})

test_that("invalid expected_environment_status errors", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$expected_environment_status[1] <- "probably"

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "expected_environment_status.*invalid values.*probably"
  )
})

test_that("missing required column errors", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$references <- NULL

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "missing required flag_taxa columns: references"
  )
})

test_that("row mismatch errors", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()[1, , drop = FALSE]

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "same number of rows"
  )
})

test_that("feature_id mismatch errors if feature_id is present", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$feature_id[2] <- "asv999"

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "feature_id.*must match"
  )
})

test_that("flag_taxa accepts phyloseq objects with taxonomy tables", {
  otu <- phyloseq::otu_table(
    matrix(
      c(1, 2),
      nrow = 1,
      dimnames = list("asv1", c("sample1", "sample2"))
    ),
    taxa_are_rows = TRUE
  )
  tax <- phyloseq::tax_table(
    matrix(
      c("Animalia", "Chordata", "Salmo salar"),
      nrow = 1,
      dimnames = list("asv1", c("kingdom", "phylum", "species"))
    )
  )
  ps <- phyloseq::phyloseq(otu, tax)

  prompt <- flag_taxa(ps, expected_environment = "marine")

  expect_type(prompt, "character")
  expect_length(prompt, 1)
  expect_true(grepl("Salmo salar", prompt, fixed = TRUE))
})

test_that("choose_query_taxon selects the most specific available rank", {
  tax_row <- c(
    kingdom = "Animalia",
    phylum = "Chordata",
    family = "Salmonidae",
    genus = "",
    species = NA_character_
  )

  query_taxon <- biohelper:::choose_query_taxon(tax_row)

  expect_equal(query_taxon$query_rank, "family")
  expect_equal(query_taxon$query_name, "Salmonidae")
})
