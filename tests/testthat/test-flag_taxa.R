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
      expected_region_status = c("known_in_region", "known_near_region"),
      recommended_action = c("retain", "review"),
      rationale = c("Example rationale 1.", "Example rationale 2."),
      references = c("Example reference 1.", "Example reference 2."),
      stringsAsFactors = FALSE
    )
  )
}

valid_taxon_evidence <- function() {
  data.frame(
    taxon_name = c("Salmo salar", "Daphnia"),
    taxon_rank = c("species", "genus"),
    source = c("literature", "literature"),
    evidence_type = c("environment", "habitat"),
    evidence_summary = c(
      "Occurs in marine and freshwater phases.",
      "Mostly freshwater, with some brackish records."
    ),
    reference = c("Example reference 1", "Example reference 2"),
    environment = c("marine; freshwater", "freshwater; brackish"),
    habitat = c("coastal water; rivers", "ponds; lakes"),
    region = c("North Atlantic", "global"),
    reference_url = c("https://example.org/salmo", "https://example.org/daphnia"),
    doi = c("10.0000/example.1", "10.0000/example.2"),
    stringsAsFactors = FALSE
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

test_that("prompt_only = TRUE asks for valid JSON only", {
  prompt <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )

  json_fields <- c(
    "feature_id",
    "query_rank",
    "query_name",
    "lineage",
    "expected_environment",
    "expected_environment_status",
    "expected_habitat",
    "expected_habitat_status",
    "expected_region",
    "expected_region_status",
    "recommended_action",
    "rationale",
    "references"
  )

  expect_true(grepl("Return valid JSON only.", prompt, fixed = TRUE))
  expect_true(grepl("JSON must be an array of objects", prompt, fixed = TRUE))
  expect_true(grepl("Do not return markdown", prompt, fixed = TRUE))
  expect_true(grepl("plain text table", prompt, fixed = TRUE))
  expect_false(grepl("Return one table only", prompt, fixed = TRUE))

  for (json_field in json_fields) {
    expect_true(grepl(paste0("- ", json_field), prompt, fixed = TRUE))
  }
})

test_that("taxon_evidence = NULL keeps current prompt behaviour", {
  prompt_default <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )
  prompt_null <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    taxon_evidence = NULL
  )

  expect_identical(prompt_null, prompt_default)
})

test_that("valid taxon_evidence is accepted", {
  prompt <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    taxon_evidence = valid_taxon_evidence()
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
})

test_that("missing required evidence columns errors clearly", {
  taxon_evidence <- valid_taxon_evidence()
  taxon_evidence$reference <- NULL

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      taxon_evidence = taxon_evidence
    ),
    "taxon_evidence.*missing required columns: reference"
  )
})

test_that("non-data.frame taxon_evidence errors clearly", {
  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      taxon_evidence = list(taxon_name = "Salmo salar")
    ),
    "`taxon_evidence` must be a data.frame"
  )
})

test_that("entirely empty required evidence columns error clearly", {
  taxon_evidence <- valid_taxon_evidence()
  taxon_evidence$evidence_summary <- ""

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      taxon_evidence = taxon_evidence
    ),
    "Required `taxon_evidence` columns.*evidence_summary"
  )
})

test_that("prompt_only = TRUE includes evidence section when taxon_evidence is provided", {
  prompt <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    taxon_evidence = valid_taxon_evidence()
  )

  expect_true(grepl("User-provided taxon evidence", prompt, fixed = TRUE))
  expect_true(grepl("Use this evidence before performing online searches.", prompt, fixed = TRUE))
  expect_true(grepl("Do not search online for taxa where the provided evidence is sufficient", prompt, fixed = TRUE))
  expect_true(grepl("Use online sources only when the provided evidence is incomplete, missing, or conflicting.", prompt, fixed = TRUE))
  expect_true(grepl("Salmo salar", prompt, fixed = TRUE))
  expect_true(grepl("Example reference 1", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE does not include evidence section when taxon_evidence is NULL", {
  prompt <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )

  expect_false(grepl("User-provided taxon evidence", prompt, fixed = TRUE))
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
      expected_region = "North Atlantic",
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
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "expected_environment_status.*invalid values.*probably"
  )
})

test_that("invalid expected_habitat_status errors", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$expected_habitat_status[1] <- "probably"

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "expected_habitat_status.*invalid values.*probably"
  )
})

test_that("invalid expected_region_status errors", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$expected_region_status[1] <- "probably"

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "expected_region_status.*invalid values.*probably"
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
      expected_region = "North Atlantic",
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
      expected_region = "North Atlantic",
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
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "feature_id.*must match"
  )
})

test_that("expected_region = NULL requires not_assessed", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "expected_region_status.*must be `not_assessed`"
  )
})

test_that("expected_region not NULL rejects not_assessed", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_result()
  result$expected_region_status[1] <- "not_assessed"

  expect_error(
    flag_taxa(
      tax,
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "expected_region_status.*must not be `not_assessed`"
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
