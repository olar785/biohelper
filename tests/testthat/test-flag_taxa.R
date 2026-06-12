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
