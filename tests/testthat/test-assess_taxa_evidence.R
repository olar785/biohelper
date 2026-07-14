assess_taxa_evidence_test_taxonomy <- function() {
  data.frame(
    feature_id = c("asv_salmo", "asv_daphnia", "asv_empty"),
    phylum = c("Chordata", "Arthropoda", NA_character_),
    class = c("Actinopteri", NA_character_, NA_character_),
    order = c(NA_character_, NA_character_, NA_character_),
    family = c("Salmonidae", NA_character_, NA_character_),
    genus = c("Salmo", "Daphnia", NA_character_),
    species = c("Salmo salar", NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
}

assess_taxa_evidence_test_evidence <- function() {
  data.frame(
    taxon_name = "Salmo salar",
    taxon_rank = "species",
    source = "literature",
    evidence_type = "environment",
    evidence_summary = "Occurs in marine and freshwater phases.",
    reference = "Example reference",
    environment = "marine; freshwater",
    habitat = "coastal water; rivers",
    region = "North Atlantic",
    reference_url = "https://example.org/salmo",
    stringsAsFactors = FALSE
  )
}

assess_taxa_evidence_rank_summary_evidence <- function() {
  data.frame(
    taxon_name = c(
      "Chordata",
      "Actinopteri",
      "Salmoniformes",
      "Salmonidae",
      "Salmo",
      "Salmo salar",
      "Salmo trutta"
    ),
    taxon_rank = c("phylum", "class", "order", "family", "genus", "species", "species"),
    source = "curated_checklist",
    evidence_type = "regional_checklist",
    evidence_summary = "Matched curated checklist evidence.",
    reference = "Curated checklist",
    environment = c("marine", "marine", "marine", "marine", "marine; brackish", "marine", "marine"),
    habitat = NA_character_,
    region = "Tropical Pacific",
    phylum = "Chordata",
    class = c(NA_character_, "Actinopteri", "Actinopteri", "Actinopteri", "Actinopteri", "Actinopteri", "Actinopteri"),
    order = c(NA_character_, NA_character_, "Salmoniformes", "Salmoniformes", "Salmoniformes", "Salmoniformes", "Salmoniformes"),
    family = c(NA_character_, NA_character_, NA_character_, "Salmonidae", "Salmonidae", "Salmonidae", "Salmonidae"),
    genus = c(NA_character_, NA_character_, NA_character_, NA_character_, "Salmo", "Salmo", "Salmo"),
    stringsAsFactors = FALSE
  )
}

assess_taxa_evidence_obis_checklist <- function(records = c(5, 20, 7, 3, 9, 11)) {
  data.frame(
    scientificName = c(
      "Chordata",
      "Actinopteri",
      "Salmoniformes",
      "Salmonidae",
      "Salmo",
      "Salmo salar"
    ),
    taxonRank = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    records = records,
    is_marine = TRUE,
    is_brackish = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    is_freshwater = FALSE,
    is_terrestrial = FALSE,
    stringsAsFactors = FALSE
  )
}

assess_taxa_evidence_local_evidence <- function(
  taxon_name,
  taxon_rank,
  environment = "marine",
  habitat = NA_character_,
  region = NA_character_,
  evidence_summary = "Curated local evidence.",
  reference = "Curated reference"
) {
  data.frame(
    taxon_name = taxon_name,
    taxon_rank = taxon_rank,
    source = "curated_checklist",
    evidence_type = "regional_checklist",
    evidence_summary = evidence_summary,
    reference = reference,
    environment = environment,
    habitat = habitat,
    region = region,
    stringsAsFactors = FALSE
  )
}

assess_taxa_evidence_ecology_evidence <- function() {
  data.frame(
    source = c("sealifebase", "fishbase", "fishbase", "fishbase"),
    spec_code = c("10", "20", NA_character_, NA_character_),
    taxon_rank = c("species", "species", "family", "subfamily"),
    taxon_name = c("Abyssogena southwardae", "Katsuwonus pelamis", "Calanidae", "Calaninae"),
    kingdom = "Metazoa",
    phylum = c("Mollusca", "Chordata", "Arthropoda", "Arthropoda"),
    class = c("Bivalvia", "Actinopterygii", "Hexanauplia", "Hexanauplia"),
    order = c(NA_character_, "Scombriformes", "Calanoida", "Calanoida"),
    family = c("Vesicomyidae", "Scombridae", "Calanidae", "Calanidae"),
    genus = c("Abyssogena", "Katsuwonus", NA_character_, NA_character_),
    species = c("southwardae", "pelamis", NA_character_, NA_character_),
    common_name = c(NA_character_, "skipjack tuna", "calanoid copepods", "calanine copepods"),
    rank_habitat_note = c(NA_character_, NA_character_, "marine planktonic copepods", "intermediate rank note"),
    rank_water_salinity = c(NA_character_, "marine", "marine", "marine"),
    rank_body_shape = c(NA_character_, "fusiform", "copepod", "copepod"),
    rank_distribution = c(NA_character_, "Pacific", "global ocean", "global ocean"),
    rank_diagnosis = c(NA_character_, "fish diagnosis", "copepod diagnosis", "intermediate diagnosis"),
    rank_remarks = c(NA_character_, "fish remarks", "copepod remarks", "intermediate remarks"),
    rank_etymology = c(NA_character_, "fish etymology", "copepod etymology", "intermediate etymology"),
    ecology_environment = c("freshwater", "marine", "marine", "marine"),
    ecology_position = c("benthic", "pelagic", "pelagic", "pelagic"),
    ecology_habitat_broad = c("benthic", "pelagic", "pelagic", "pelagic"),
    ecology_depth_zone = c("abyssal_possible", "shallow_only", "shallow_only", "shallow_only"),
    depth_min = c(2400, 0, 0, 0),
    depth_max = c(4200, 260, 600, 600),
    common_depth_min = c(2500, NA_real_, NA_real_, NA_real_),
    common_depth_max = c(4000, NA_real_, NA_real_, NA_real_),
    n_species = c(NA_real_, NA_real_, 42, 12),
    dominant_habitat_broad = c(NA_character_, NA_character_, "pelagic", "pelagic"),
    dominant_habitat_broad_prop = c(NA_real_, NA_real_, 0.95, 0.9),
    dominant_depth_zone = c(NA_character_, NA_character_, "shallow_only", "shallow_only"),
    dominant_depth_zone_prop = c(NA_real_, NA_real_, 0.8, 0.8),
    ecology_evidence_strength = c(NA_character_, NA_character_, "high", "high"),
    stringsAsFactors = FALSE
  )
}

assess_taxa_evidence_details <- function(...) {
  assess_taxa_evidence(...)$details
}

test_that("assess_taxa_evidence is exported with a deterministic-only interface", {
  expect_true("assess_taxa_evidence" %in% getNamespaceExports("biohelper"))

  forbidden_args <- c(
    "chat",
    "allow_llm_tools",
    "judgement_mode",
    "prompt_only",
    "prompt_path",
    "prompt_tools",
    "tool_requirement",
    "llm_result",
    "mock_llm_result",
    "max_tool_taxa",
    "tool_batch_size",
    "max_evidence_rows_per_query",
    "ccz_evidence_path"
  )
  expect_false(any(forbidden_args %in% names(formals(assess_taxa_evidence))))
  expect_true("taxon_evidence_path" %in% names(formals(assess_taxa_evidence)))
  expect_true("ecology_evidence" %in% names(formals(assess_taxa_evidence)))
  expect_true("ecology_evidence_path" %in% names(formals(assess_taxa_evidence)))
  expect_true("expected_habitat_type" %in% names(formals(assess_taxa_evidence)))
  expect_false("output" %in% names(formals(assess_taxa_evidence)))
})

test_that("assess_taxa_evidence runs without chat and returns deterministic columns", {
  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = assess_taxa_evidence_test_evidence()
  )

  expected_columns <- c(
    "feature_id",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale",
    "references",
    "evidence_basis",
    "evidence_sources",
    "evidence_summary"
  )
  expect_true(all(expected_columns %in% colnames(output)))
  expect_false(any(grepl("^llm_", colnames(output))))
  expect_false(any(c(
    "query_rank",
    "query_name",
    "matched_rank",
    "matched_name"
  ) %in% colnames(output)))
  expect_equal(output$feature_id, c("asv_salmo", "asv_daphnia", "asv_empty"))
  expect_equal(output$taxon_name[[1]], "Salmo salar")
  expect_equal(output$taxon_rank[[1]], "species")
  expect_equal(output$recommended_action[[1]], "retain")
  expect_equal(output$evidence_sources[[1]], "literature")
})

test_that("assess_taxa_evidence matches flag_taxa evidence-only deterministic output", {
  tax <- assess_taxa_evidence_test_taxonomy()
  evidence <- assess_taxa_evidence_test_evidence()

  assessed <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )
  flag_output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  comparable_columns <- c(
    "feature_id",
    "taxon_name",
    "taxon_rank",
    "lineage",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale",
    "references",
    "evidence_basis",
    "evidence_sources",
    "evidence_summary"
  )
  assessed_comparable <- assessed[1:2, comparable_columns, drop = FALSE]
  flag_comparable <- flag_output[1:2, comparable_columns, drop = FALSE]
  expect_equal(assessed_comparable, flag_comparable)
})

test_that("assess_taxa_evidence treats missing evidence as uncertainty, not exclusion", {
  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1:2, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "deep sea"
  )

  expect_equal(output$recommended_action, rep("flag_for_review", nrow(output)))
  expect_false(any(output$recommended_action == "exclude"))
  expect_true(all(output$expected_environment_status == "unknown"))
  # An expected habitat was supplied, but no deterministic habitat evidence was matched.
  expect_true(all(output$expected_habitat_status == "no_known_habitat_evidence"))
  expect_true(all(grepl("Missing evidence", output$rationale, fixed = TRUE)))
  expect_equal(output$evidence_sources, rep("none", nrow(output)))
})

test_that("summarise_obis_checklist summarises mixed ranks and environments", {
  chk <- rbind(
    assess_taxa_evidence_obis_checklist(),
    data.frame(
      scientificName = "Salmo trutta",
      taxonRank = "Species",
      records = 30,
      is_marine = TRUE,
      is_brackish = TRUE,
      is_freshwater = FALSE,
      is_terrestrial = FALSE,
      stringsAsFactors = FALSE
    )
  )

  summary <- biohelper:::summarise_obis_checklist(
    chk,
    query_name = "Salmo salar",
    region = "North Atlantic"
  )

  expect_equal(summary$obis_region_status, "obis_records_found")
  expect_equal(summary$obis_taxa_returned, 7)
  expect_equal(summary$obis_total_records_including_descendants, sum(chk$records))
  expect_equal(summary$obis_environment, "marine; brackish")
  expect_equal(
    colnames(summary)[match("obis_n_phylum", colnames(summary)):match("obis_taxa_species", colnames(summary))],
    c(
      "obis_n_phylum", "obis_taxa_phylum",
      "obis_n_class", "obis_taxa_class",
      "obis_n_order", "obis_taxa_order",
      "obis_n_family", "obis_taxa_family",
      "obis_n_genus", "obis_taxa_genus",
      "obis_n_species", "obis_taxa_species"
    )
  )
  expect_equal(summary$obis_n_species, 2)
  expect_equal(summary$obis_taxa_species, "Salmo trutta | Salmo salar")
})

test_that("summarise_obis_checklist returns empty summary for no records", {
  summary <- biohelper:::summarise_obis_checklist(
    data.frame(),
    query_name = "Daphnia",
    region = "North Atlantic"
  )

  expect_equal(summary$obis_region_status, "no_obis_records")
  expect_equal(summary$obis_taxa_returned, 0)
  expect_equal(summary$obis_total_records_including_descendants, 0)
  expect_true(is.na(summary$obis_environment))
  rank_count_columns <- paste0("obis_n_", c("phylum", "class", "order", "family", "genus", "species"))
  rank_taxa_columns <- paste0("obis_taxa_", c("phylum", "class", "order", "family", "genus", "species"))
  expect_true(all(unlist(summary[rank_count_columns]) == 0))
  expect_true(all(is.na(unlist(summary[rank_taxa_columns]))))
})

test_that("assess_taxa_evidence uses direct OBIS checklist records as known regional evidence", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) {
      expect_equal(scientificname, "Salmo salar")
      expect_equal(geometry, "POLYGON ((0 0, 1 0, 1 1, 0 0))")
      assess_taxa_evidence_obis_checklist()
    }
  )

  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "North Atlantic",
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))",
    obis_region_relation = "direct"
  )

  expect_equal(output$expected_region_status, "known_in_region")
  expect_equal(output$obis_region_status, "obis_records_found")
  expect_equal(output$obis_taxa_returned, 6)
  expect_equal(output$obis_environment, "marine; brackish")
  expect_true(grepl("OBIS checklist query for Salmo salar in North Atlantic", output$evidence_summary, fixed = TRUE))
  expect_true(grepl("obis", output$evidence_sources, fixed = TRUE))
  expect_false(output$recommended_action == "exclude")
})

test_that("assess_taxa_evidence uses supporting OBIS checklist records as plausible regional evidence", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) assess_taxa_evidence_obis_checklist()
  )

  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_region = "Cook Islands",
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))",
    obis_region_label = "Tropical Pacific",
    obis_region_relation = "supporting"
  )

  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_equal(output$obis_region_status, "obis_records_found")
  expect_false(output$recommended_action == "exclude")
})

test_that("assess_taxa_evidence applies direct and supporting user evidence region relations", {
  tax <- assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE]
  evidence <- assess_taxa_evidence_test_evidence()

  direct <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "Cook Islands",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )
  supporting <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "Cook Islands",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "supporting"
  )

  expect_equal(direct$expected_region_status, "known_in_region")
  expect_equal(supporting$expected_region_status, "plausible_in_region")
  expect_equal(direct$taxon_evidence_region, "North Atlantic")
  expect_equal(supporting$taxon_evidence_region_relation, "supporting")
})

test_that("assess_taxa_evidence accepts taxon_evidence_path", {
  evidence_path <- tempfile(fileext = ".csv")
  utils::write.csv(
    assess_taxa_evidence_test_evidence(),
    evidence_path,
    row.names = FALSE
  )

  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence_path = evidence_path,
    taxon_evidence_region_relation = "direct"
  )

  expect_equal(output$expected_region_status, "known_in_region")
  expect_equal(output$taxon_evidence_taxa_returned, 1)
  expect_equal(output$taxon_evidence_taxa_species, "Salmo salar")
})

test_that("assess_taxa_evidence returns user evidence summary columns in rank order", {
  tax <- data.frame(
    feature_id = "asv_chordata",
    phylum = "Chordata",
    stringsAsFactors = FALSE
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_region = "Cook Islands",
    taxon_evidence = assess_taxa_evidence_rank_summary_evidence(),
    taxon_evidence_region_label = "Tropical Pacific",
    taxon_evidence_region_relation = "supporting"
  )

  expected_summary_columns <- c(
    "taxon_evidence_region",
    "taxon_evidence_region_relation",
    "taxon_evidence_taxa_returned",
    "taxon_evidence_environment",
    "taxon_evidence_n_phylum",
    "taxon_evidence_taxa_phylum",
    "taxon_evidence_n_class",
    "taxon_evidence_taxa_class",
    "taxon_evidence_n_order",
    "taxon_evidence_taxa_order",
    "taxon_evidence_n_family",
    "taxon_evidence_taxa_family",
    "taxon_evidence_n_genus",
    "taxon_evidence_taxa_genus",
    "taxon_evidence_n_species",
    "taxon_evidence_taxa_species"
  )
  expect_true(all(expected_summary_columns %in% colnames(output)))
  expect_equal(
    colnames(output)[match("taxon_evidence_region", colnames(output)):match("taxon_evidence_taxa_species", colnames(output))],
    expected_summary_columns
  )
  expect_equal(output$taxon_evidence_region, "Tropical Pacific")
  expect_equal(output$taxon_evidence_region_relation, "supporting")
  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_equal(output$taxon_evidence_taxa_returned, 7)
  expect_equal(output$taxon_evidence_environment, "marine; brackish")
  expect_equal(output$taxon_evidence_n_phylum, 1)
  expect_equal(output$taxon_evidence_taxa_phylum, "Chordata")
  expect_equal(output$taxon_evidence_n_species, 2)
  expect_equal(output$taxon_evidence_taxa_species, "Salmo salar | Salmo trutta")
})

test_that("assess_taxa_evidence fills user evidence summary defaults when no rows match", {
  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[2, , drop = FALSE],
    expected_environment = "marine",
    expected_region = "North Atlantic",
    taxon_evidence = assess_taxa_evidence_test_evidence()
  )

  expect_equal(output$taxon_evidence_taxa_returned, 0)
  expect_true(is.na(output$taxon_evidence_environment))
  expect_equal(output$taxon_evidence_n_species, 0)
  expect_true(is.na(output$taxon_evidence_taxa_species))
})

test_that("assess_taxa_evidence output begins with primary interpretation columns", {
  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = assess_taxa_evidence_test_evidence(),
    taxon_evidence_region_relation = "direct"
  )

  expected_prefix <- c(
    "feature_id",
    "taxon_rank",
    "taxon_name",
    "lineage",
    "worms_environment",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale",
    "references",
    "evidence_basis",
    "evidence_sources",
    "evidence_summary"
  )
  expect_equal(colnames(output)[seq_along(expected_prefix)], expected_prefix)
  expect_gt(match("taxon_evidence_region", colnames(output)), match("evidence_summary", colnames(output)))
  expect_gt(match("expected_environment", colnames(output)), match("obis_taxa_species", colnames(output), nomatch = 0))
})

test_that("OBIS absence does not override supporting user evidence or trigger exclusion", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) data.frame()
  )

  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "Cook Islands",
    taxon_evidence = assess_taxa_evidence_test_evidence(),
    taxon_evidence_region_relation = "supporting",
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))"
  )

  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_equal(output$obis_region_status, "no_obis_records")
  expect_false(output$recommended_action == "exclude")
})

test_that("OBIS absence alone is no distribution evidence, not exclusion", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) data.frame()
  )

  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "North Atlantic",
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))"
  )

  expect_equal(output$expected_region_status, "no_distribution_evidence")
  expect_equal(output$obis_region_status, "no_obis_records")
  expect_false(output$recommended_action == "exclude")
  expect_false(output$expected_environment_status == "incompatible")
})

test_that("OBIS requested without geometry warns and skips OBIS", {
  output <- NULL
  expect_warning(
    output <- assess_taxa_evidence_details(
      assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
      expected_environment = "marine",
      expected_region = "North Atlantic",
      evidence_sources = "obis"
    ),
    "OBIS regional evidence requires `obis_geometry`"
  )

  expect_false("obis_region_status" %in% colnames(output))
  expect_equal(output$expected_region_status, "no_distribution_evidence")
})

test_that("OBIS and taxon evidence region relations validate direct or supporting only", {
  expect_error(
    assess_taxa_evidence(
      assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
      expected_environment = "marine",
      obis_region_relation = "unknown"
    ),
    "'arg' should be one of"
  )
  expect_error(
    assess_taxa_evidence(
      assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
      expected_environment = "marine",
      taxon_evidence_region_relation = "unknown"
    ),
    "'arg' should be one of"
  )
})

test_that("compatible environment with no habitat evidence and plausible region remains review", {
  tax <- data.frame(
    feature_id = "asv_delphinus",
    phylum = "Chordata",
    class = "Mammalia",
    order = "Cetacea",
    family = "Delphinidae",
    genus = "Delphinus",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Delphinus",
    taxon_rank = "genus",
    environment = "marine",
    habitat = NA_character_,
    region = "Tropical Pacific"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_label = "Tropical Pacific",
    taxon_evidence_region_relation = "supporting"
  )

  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_equal(output$ecological_status, "plausible")
  expect_equal(output$occurrence_interpretation, "plausible_resident")
  # Region evidence supports plausibility, but expected-habitat evidence is still missing.
  expect_equal(output$recommended_action, "flag_for_review")
  expect_true(grepl("No direct expected-habitat evidence was matched", output$rationale, fixed = TRUE))
  expect_true(grepl("Region status: plausible_in_region", output$rationale, fixed = TRUE))
})

test_that("compatible environment with no habitat evidence and known region remains review", {
  tax <- data.frame(
    feature_id = "asv_delphinus",
    phylum = "Chordata",
    class = "Mammalia",
    genus = "Delphinus",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Delphinus",
    taxon_rank = "genus",
    environment = "marine",
    habitat = NA_character_,
    region = "North Atlantic"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "North Atlantic",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )

  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_equal(output$expected_region_status, "known_in_region")
  # Missing expected-habitat evidence is uncertainty and should not be promoted to retain.
  expect_equal(output$recommended_action, "flag_for_review")
  expect_equal(output$occurrence_interpretation, "plausible_resident")
})

test_that("broad compatible taxa with no habitat evidence remain review when region is plausible", {
  tax <- data.frame(
    feature_id = "asv_calanoida",
    phylum = "Arthropoda",
    class = "Copepoda",
    order = "Calanoida",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Calanoida",
    taxon_rank = "order",
    environment = "marine; freshwater",
    habitat = NA_character_,
    region = "Tropical Pacific",
    evidence_summary = "Calanoida includes marine and freshwater members."
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_label = "Tropical Pacific",
    taxon_evidence_region_relation = "supporting"
  )

  expect_equal(output$expected_environment_status, "mixed_within_rank")
  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_equal(output$recommended_action, "flag_for_review")
  expect_true(grepl("missing habitat or region evidence is treated cautiously", output$rationale, fixed = TRUE))
})

test_that("supporting OBIS region evidence does not retain without habitat evidence", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) {
      expect_equal(scientificname, "Delphinus")
      assess_taxa_evidence_obis_checklist(records = c(2, 4, 3, 2, 5, 0))
    }
  )
  tax <- data.frame(
    feature_id = "asv_delphinus",
    phylum = "Chordata",
    class = "Mammalia",
    genus = "Delphinus",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Delphinus",
    taxon_rank = "genus",
    environment = "marine",
    habitat = NA_character_
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "supporting",
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))",
    obis_region_label = "Tropical Pacific",
    obis_region_relation = "supporting"
  )

  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_equal(output$obis_region_status, "obis_records_found")
  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_equal(output$recommended_action, "flag_for_review")
})

test_that("photosynthetic deep-sea mismatch remains possible exclusion", {
  tax <- data.frame(
    feature_id = "asv_florideophyceae",
    phylum = "Rhodophyta",
    class = "Florideophyceae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Florideophyceae",
    taxon_rank = "class",
    environment = "marine",
    habitat = NA_character_,
    region = "Clarion-Clipperton Zone",
    evidence_summary = "Florideophyceae are photosynthetic red algal taxa."
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )

  expect_equal(output$expected_habitat_status, "incompatible")
  expect_equal(output$recommended_action, "flag_possible_exclusion")
  expect_equal(output$occurrence_interpretation, "possible_transient_or_allochthonous")
})

test_that("missing taxonomy remains review with no deprecated habitat status", {
  tax <- data.frame(
    feature_id = "asv_empty",
    kingdom = "Eukaryota",
    phylum = NA_character_,
    class = NA_character_,
    order = NA_character_,
    family = NA_character_,
    genus = NA_character_,
    species = NA_character_,
    stringsAsFactors = FALSE
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone"
  )

  expect_equal(output$recommended_action, "flag_for_review")
  expect_true(output$expected_habitat_status %in% c("unknown", "insufficient_taxonomic_resolution"))
  expect_false(any(output$expected_habitat_status == "no_known_habitat_evidence"))
})

test_that("returned habitat statuses use explicit no-known evidence wording", {
  output <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1:2, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "deep sea"
  )

  # The expected habitat was assessed, and no matching habitat evidence was found.
  expect_true(all(output$expected_habitat_status == "no_known_habitat_evidence"))
  expect_false(any(output$expected_habitat_status == "unknown"))
})

test_that("terrestrial-only taxa in expected marine context are excluded", {
  tax <- data.frame(
    feature_id = "asv_bovidae",
    phylum = "Chordata",
    class = "Mammalia",
    order = "Artiodactyla",
    family = "Bovidae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Bovidae",
    taxon_rank = "family",
    environment = "terrestrial",
    habitat = NA_character_,
    region = NA_character_,
    evidence_summary = "Bovidae are terrestrial mammals."
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence
  )

  expect_equal(output$expected_environment_status, "incompatible")
  expect_equal(output$ecological_status, "incompatible_resident")
  expect_equal(output$occurrence_interpretation, "likely_contaminant_or_misassignment")
  expect_equal(output$recommended_action, "exclude")
  expect_true(grepl("terrestrial", output$rationale, fixed = TRUE))
  expect_true(grepl("positively contradicts the expected environment", output$rationale, fixed = TRUE))
})

test_that("freshwater-only taxa in expected marine context are excluded", {
  tax <- data.frame(
    feature_id = "asv_daphnia",
    phylum = "Arthropoda",
    class = "Branchiopoda",
    genus = "Daphnia",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Daphnia",
    taxon_rank = "genus",
    environment = "freshwater",
    habitat = NA_character_,
    region = NA_character_,
    evidence_summary = "Daphnia evidence is freshwater-only."
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = evidence
  )

  expect_equal(output$expected_environment_status, "incompatible")
  expect_equal(output$recommended_action, "exclude")
  expect_false(output$expected_environment_status == "mixed_within_rank")
})

test_that("mixed environment including marine remains mixed and is not excluded by environment alone", {
  tax <- data.frame(
    feature_id = "asv_calanoida",
    phylum = "Arthropoda",
    class = "Copepoda",
    order = "Calanoida",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Calanoida",
    taxon_rank = "order",
    environment = "marine; brackish; freshwater",
    habitat = NA_character_,
    region = "Tropical Pacific"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_label = "Tropical Pacific",
    taxon_evidence_region_relation = "supporting"
  )

  expect_equal(output$expected_environment_status, "mixed_within_rank")
  expect_equal(output$recommended_action, "flag_for_review")
  expect_false(output$recommended_action == "exclude")
})

test_that("OBIS records do not override trusted terrestrial environment evidence", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) {
      expect_equal(scientificname, "Bovidae")
      data.frame(
        scientificName = "Bovidae",
        taxonRank = "Family",
        records = 4,
        is_marine = TRUE,
        is_brackish = FALSE,
        is_freshwater = FALSE,
        is_terrestrial = TRUE,
        stringsAsFactors = FALSE
      )
    }
  )
  tax <- data.frame(
    feature_id = "asv_bovidae",
    phylum = "Chordata",
    class = "Mammalia",
    family = "Bovidae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Bovidae",
    taxon_rank = "family",
    environment = "terrestrial",
    habitat = NA_character_,
    region = NA_character_
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))",
    obis_region_label = "Tropical Pacific",
    obis_region_relation = "supporting"
  )

  expect_equal(output$obis_environment, "marine; terrestrial")
  expect_equal(output$expected_environment_status, "incompatible")
  expect_equal(output$recommended_action, "exclude")
})

test_that("OBIS environment flags are fallback only when trusted environment evidence is missing", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) {
      expect_equal(scientificname, "Delphinus")
      data.frame(
        scientificName = "Delphinus",
        taxonRank = "Genus",
        records = 8,
        is_marine = TRUE,
        is_brackish = FALSE,
        is_freshwater = FALSE,
        is_terrestrial = FALSE,
        stringsAsFactors = FALSE
      )
    }
  )
  tax <- data.frame(
    feature_id = "asv_delphinus",
    phylum = "Chordata",
    class = "Mammalia",
    genus = "Delphinus",
    stringsAsFactors = FALSE
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    evidence_sources = "obis",
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))",
    obis_region_label = "Tropical Pacific",
    obis_region_relation = "supporting"
  )

  expect_equal(output$expected_environment_status, "compatible")
  expect_true(grepl("OBIS environment fallback", output$evidence_summary, fixed = TRUE))
  expect_true(grepl("OBIS environment flags", output$rationale, fixed = TRUE))
})

test_that("marine environment with unknown habitat and no region remains review not exclusion", {
  tax <- data.frame(
    feature_id = "asv_abyssogena",
    phylum = "Mollusca",
    class = "Bivalvia",
    genus = "Abyssogena",
    species = "Abyssogena southwardae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Abyssogena southwardae",
    taxon_rank = "species",
    environment = "marine",
    habitat = NA_character_,
    region = NA_character_
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence
  )

  expect_equal(output$expected_environment_status, "compatible")
  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_equal(output$expected_region_status, "no_distribution_evidence")
  expect_equal(output$recommended_action, "flag_for_review")
  expect_false(output$recommended_action == "exclude")
})

test_that("final rationales match merged statuses and actions", {
  tax <- data.frame(
    feature_id = "asv_calanoida",
    phylum = "Arthropoda",
    class = "Copepoda",
    order = "Calanoida",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Calanoida",
    taxon_rank = "order",
    environment = "marine; freshwater",
    habitat = NA_character_,
    region = "Tropical Pacific"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_label = "Tropical Pacific",
    taxon_evidence_region_relation = "supporting"
  )

  expect_equal(output$recommended_action, "flag_for_review")
  expect_equal(output$expected_region_status, "plausible_in_region")
  expect_true(grepl("missing habitat or region evidence is treated cautiously", output$rationale, fixed = TRUE))
  expect_false(grepl("Region status: known_elsewhere_only", output$rationale, fixed = TRUE))
  expect_true(grepl("Region status: plausible_in_region", output$rationale, fixed = TRUE))
  expect_true(grepl("Habitat status: no_known_habitat_evidence", output$rationale, fixed = TRUE))
  expect_true(grepl("no_known_habitat_evidence", output$rationale, fixed = TRUE))
})

test_that("FishBase/SeaLifeBase ecology evidence refines habitat but not environment or region", {
  tax <- data.frame(
    feature_id = "asv_abyssogena",
    phylum = "Mollusca",
    class = "Bivalvia",
    family = "Vesicomyidae",
    genus = "Abyssogena",
    species = "Abyssogena southwardae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Abyssogena southwardae",
    taxon_rank = "species",
    environment = "marine",
    habitat = NA_character_,
    region = NA_character_
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    ecology_evidence = assess_taxa_evidence_ecology_evidence()
  )

  expect_equal(output$expected_environment_status, "compatible")
  expect_equal(output$expected_region_status, "no_distribution_evidence")
  expect_equal(output$expected_habitat_status, "possible_likely")
  expect_equal(output$ecology_evidence_sources, "sealifebase")
  expect_equal(output$ecology_evidence_depth_zone, "abyssal_possible")
  expect_true(grepl("FishBase/SeaLifeBase ecology evidence", output$evidence_summary, fixed = TRUE))
  expect_true(grepl("does not modify expected_environment_status or expected_region_status", output$rationale, fixed = TRUE))
})

test_that("FishBase/SeaLifeBase ecology evidence can support higher-rank habitat by lineage", {
  tax <- data.frame(
    feature_id = "asv_scombridae",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Scombriformes",
    family = "Scombridae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Scombridae",
    taxon_rank = "family",
    environment = "marine",
    habitat = NA_character_,
    region = "Tropical Pacific"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "pelagic water column",
    expected_region = "Tropical Pacific",
    taxon_evidence = evidence,
    ecology_evidence = assess_taxa_evidence_ecology_evidence()
  )

  expect_equal(output$expected_environment_status, "compatible")
  expect_equal(output$expected_region_status, "known_in_region")
  expect_equal(output$expected_habitat_status, "compatible")
  expect_equal(output$recommended_action, "retain")
  expect_true(grepl("Katsuwonus pelamis", output$ecology_evidence_taxa, fixed = TRUE))
  expect_equal(output$ecology_evidence_habitat_broad, "pelagic")
})

test_that("FishBase/SeaLifeBase ecology mismatch flags habitat without changing environment or region", {
  tax <- data.frame(
    feature_id = "asv_calanidae",
    phylum = "Arthropoda",
    class = "Hexanauplia",
    order = "Calanoida",
    family = "Calanidae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Calanidae",
    taxon_rank = "family",
    environment = "marine",
    habitat = NA_character_,
    region = "Clarion-Clipperton Zone"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    ecology_evidence = assess_taxa_evidence_ecology_evidence()
  )

  expect_equal(output$expected_environment_status, "compatible")
  expect_equal(output$expected_region_status, "known_in_region")
  expect_equal(output$expected_habitat_status, "possible_unlikely")
  expect_equal(output$ecological_status, "unlikely_resident")
  expect_equal(output$recommended_action, "flag_possible_exclusion")
  expect_false(grepl("Calaninae", output$ecology_evidence_taxa, fixed = TRUE))
})

test_that("FishBase/SeaLifeBase ecology evidence accepts builder list and path inputs", {
  ecology <- assess_taxa_evidence_ecology_evidence()
  ecology <- ecology[ecology$taxon_rank != "subfamily", , drop = FALSE]
  db <- list(combined = ecology)
  path <- tempfile(fileext = ".rds")
  saveRDS(db, path)

  tax <- data.frame(
    feature_id = "asv_scombridae",
    phylum = "Chordata",
    class = "Actinopterygii",
    family = "Scombridae",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Scombridae",
    taxon_rank = "family",
    environment = "marine",
    region = "Tropical Pacific"
  )

  from_list <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "pelagic water column",
    expected_region = "Tropical Pacific",
    taxon_evidence = evidence,
    ecology_evidence = db
  )
  from_path <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "pelagic water column",
    expected_region = "Tropical Pacific",
    taxon_evidence = evidence,
    ecology_evidence_path = path
  )

  expect_equal(from_list$ecology_evidence_taxa, from_path$ecology_evidence_taxa)
  expect_equal(from_path$expected_environment_status, "compatible")
  expect_equal(from_path$expected_region_status, "known_in_region")
})

test_that("assess_taxa_evidence always returns a list with compact summary", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) assess_taxa_evidence_obis_checklist()
  )
  output <- assess_taxa_evidence(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = assess_taxa_evidence_test_evidence(),
    ecology_evidence = assess_taxa_evidence_ecology_evidence(),
    evidence_sources = c("obis", "ecology"),
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))"
  )
  summary <- output$summary

  expect_type(output, "list")
  expect_equal(names(output), c("summary", "details", "worms", "obis", "ecology"))
  expect_true(all(c(
    "feature_id",
    "taxon_name",
    "ecology_evidence_common_name",
    "taxon_rank",
    "lineage",
    "worms_environment",
    "ecology_evidence_habitat_broad",
    "ecology_evidence_depth_zone",
    "ecology_evidence_depth_min",
    "ecology_evidence_depth_max",
    "obis_region_status",
    "obis_taxa_returned",
    "obis_total_records_including_descendants",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "recommended_action",
    "rationale",
    "references"
  ) %in% colnames(summary)))
  expect_equal(match("lineage", colnames(summary)), match("taxon_rank", colnames(summary)) + 1L)
  expect_false("query_id" %in% colnames(summary))
  expect_false(any(c("obis_n_rank", "obis_taxa_rank", "obis_n_phylum", "obis_taxa_species") %in% colnames(summary)))
  expect_false(any(c("taxon_evidence_n_phylum", "taxon_evidence_taxa_species") %in% colnames(summary)))
  expect_false(any(c(
    "ecology_evidence_dominant_habitat_broad",
    "ecology_evidence_dominant_habitat_broad_prop",
    "ecology_evidence_dominant_depth_zone",
    "ecology_evidence_dominant_depth_zone_prop"
  ) %in% colnames(summary)))
  expect_false("ecology_evidence_position" %in% colnames(summary))
  expect_true("ecology_evidence_habitat_broad" %in% colnames(summary))
  expect_lt(nchar(summary$rationale), 120)
})

test_that("returned list exposes cleaned details and source-specific evidence tables", {
  testthat::local_mocked_bindings(
    .require_robis = function() TRUE,
    .obis_checklist = function(scientificname, geometry) assess_taxa_evidence_obis_checklist()
  )
  details <- assess_taxa_evidence_details(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = assess_taxa_evidence_test_evidence(),
    ecology_evidence = assess_taxa_evidence_ecology_evidence(),
    evidence_sources = c("obis", "ecology"),
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))"
  )
  listed <- assess_taxa_evidence(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic",
    taxon_evidence = assess_taxa_evidence_test_evidence(),
    ecology_evidence = assess_taxa_evidence_ecology_evidence(),
    evidence_sources = c("obis", "ecology"),
    obis_geometry = "POLYGON ((0 0, 1 0, 1 1, 0 0))"
  )

  expect_gt(ncol(details), ncol(listed$summary))
  expect_false("query_id" %in% colnames(details))
  expect_false("ecology_evidence_position" %in% colnames(details))
  expect_true(all(c("summary", "details", "worms", "obis", "ecology") %in% names(listed)))
  expect_true("obis_n_species" %in% colnames(listed$obis))
  expect_true("obis_n_rank" %in% colnames(listed$obis))
  expect_true("ecology_evidence_dominant_depth_zone" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_dominant_depth_zone_prop" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_rank_water_salinity" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_rank_body_shape" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_rank_distribution" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_rank_diagnosis" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_rank_remarks" %in% colnames(listed$ecology))
  expect_true("ecology_evidence_rank_etymology" %in% colnames(listed$ecology))
})

test_that("ecology evidence summaries de-duplicate terms and use single dominant values", {
  tax <- data.frame(
    feature_id = "asv_testus",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Testiformes",
    family = "Testidae",
    genus = "Testus",
    species = "Testus alpha",
    stringsAsFactors = FALSE
  )
  ecology <- data.frame(
    source = "fishbase",
    taxon_rank = "species",
    taxon_name = c("Testus alpha", "Testus alpha"),
    ecology_environment = c("marine", "marine"),
    ecology_habitat_broad = c(
      "benthic; host_associated; pelagic; benthic",
      "host_associated; unknown; pelagic"
    ),
    ecology_depth_zone = c(
      "bathyal_possible; shallow_only",
      "bathyal_possible; abyssal_possible"
    ),
    stringsAsFactors = FALSE
  )

  output <- assess_taxa_evidence(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    ecology_evidence = ecology
  )
  ecology_out <- output$ecology

  expect_equal(
    ecology_out$ecology_evidence_habitat_broad,
    "benthic; host_associated; pelagic; unknown"
  )
  expect_false(grepl(";", ecology_out$ecology_evidence_dominant_habitat_broad, fixed = TRUE))
  expect_false(grepl(";", ecology_out$ecology_evidence_dominant_depth_zone, fixed = TRUE))
  expect_equal(ecology_out$ecology_evidence_dominant_habitat_broad, "host_associated")
  expect_equal(ecology_out$ecology_evidence_dominant_habitat_broad_prop, 2 / 6)
  expect_equal(ecology_out$ecology_evidence_dominant_depth_zone, "bathyal_possible")
  expect_equal(ecology_out$ecology_evidence_dominant_depth_zone_prop, 2 / 4)
})

test_that("dominant ecology helper is deterministic and handles all missing values", {
  tied <- biohelper:::.assess_dominant_semicolon_value(c("pelagic", "benthic"))
  missing <- biohelper:::.assess_dominant_semicolon_value(c(NA_character_, ""))

  expect_equal(tied$value, "benthic")
  expect_equal(tied$prop, 0.5)
  expect_true(is.na(missing$value))
  expect_true(is.na(missing$prop))
})

test_that("summary common name is a single best matched value", {
  tax <- data.frame(
    feature_id = "asv_skipjack",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Scombriformes",
    family = "Scombridae",
    genus = "Katsuwonus",
    species = "Katsuwonus pelamis",
    stringsAsFactors = FALSE
  )
  ecology <- data.frame(
    source = "fishbase",
    taxon_rank = c("family", "genus", "species"),
    taxon_name = c("Scombridae", "Katsuwonus", "Katsuwonus pelamis"),
    family = "Scombridae",
    genus = "Katsuwonus",
    species = c(NA_character_, NA_character_, "pelamis"),
    common_name = c("mackerels", "skipjacks", "skipjack tuna"),
    ecology_environment = "marine",
    ecology_habitat_broad = "pelagic",
    ecology_depth_zone = "shallow_only",
    stringsAsFactors = FALSE
  )

  output <- assess_taxa_evidence(
    tax,
    expected_environment = "marine",
    ecology_evidence = ecology
  )

  expect_equal(output$summary$ecology_evidence_common_name, "skipjack tuna")
  expect_false(grepl(";", output$summary$ecology_evidence_common_name, fixed = TRUE))
})

test_that("higher-rank ecology metadata is selected only from exact native rank rows", {
  tax <- data.frame(
    feature_id = "asv_carcharhinus",
    phylum = "Chordata",
    class = "Chondrichthyes",
    order = "Carcharhiniformes",
    family = "Carcharhinidae",
    genus = "Carcharhinus",
    species = NA_character_,
    stringsAsFactors = FALSE
  )
  ecology <- data.frame(
    source = "fishbase",
    taxon_rank = c("genus", "species"),
    taxon_name = c("Carcharhinus", "Carcharhinus acronotus"),
    family = "Carcharhinidae",
    genus = "Carcharhinus",
    species = c(NA_character_, "acronotus"),
    common_name = c(NA_character_, "Blacknose shark"),
    rank_water_salinity = c("marine", NA_character_),
    rank_body_shape = c(NA_character_, "species shape should not leak"),
    rank_habitat_note = c("native genus habitat note", NA_character_),
    rank_etymology = c("native genus etymology", NA_character_),
    ecology_environment = "marine",
    ecology_habitat_broad = "pelagic",
    ecology_depth_zone = "shallow_only",
    stringsAsFactors = FALSE
  )

  output <- assess_taxa_evidence(
    tax,
    expected_environment = "marine",
    ecology_evidence = ecology
  )

  expect_true(is.na(output$summary$ecology_evidence_common_name))
  expect_equal(output$ecology$ecology_evidence_rank_water_salinity, "marine")
  expect_true(is.na(output$ecology$ecology_evidence_rank_body_shape))
  expect_equal(output$ecology$ecology_evidence_rank_habitat_note, "native genus habitat note")
  expect_equal(output$ecology$ecology_evidence_rank_etymology, "native genus etymology")
})

test_that("expected habitat type uses the controlled vocabulary and inference", {
  expect_equal(
    biohelper:::.validate_assess_expected_habitat_type(NULL, "deep-sea"),
    "deep_sea_benthic"
  )
  expect_equal(
    biohelper:::.validate_assess_expected_habitat_type("deep_sea_pelagic", NULL),
    "deep_sea_pelagic"
  )
  expect_equal(
    biohelper:::.validate_assess_expected_habitat_type(NULL, NULL),
    "unspecified"
  )
  expect_error(
    biohelper:::.validate_assess_expected_habitat_type("deep-sea", NULL),
    "`expected_habitat_type` must be one of"
  )
})

test_that("deep-sea benthic Aves and Reptilia are possible allochthonous but mammals are not blanket-flagged", {
  taxa <- data.frame(
    feature_id = c("asv_bird", "asv_reptile", "asv_mammal"),
    phylum = "Chordata",
    class = c("Aves", "Reptilia", "Mammalia"),
    order = c("Charadriiformes", "Testudines", "Cetacea"),
    family = c("Laridae", "Cheloniidae", "Ziphiidae"),
    genus = c("Larus", "Chelonia", "Ziphius"),
    species = c("Larus dominicanus", "Chelonia mydas", "Ziphius cavirostris"),
    stringsAsFactors = FALSE
  )
  evidence <- dplyr::bind_rows(
    assess_taxa_evidence_local_evidence("Larus dominicanus", "species", environment = "marine", region = "Clarion-Clipperton Zone"),
    assess_taxa_evidence_local_evidence("Chelonia mydas", "species", environment = "marine", region = "Clarion-Clipperton Zone"),
    assess_taxa_evidence_local_evidence("Ziphius cavirostris", "species", environment = "marine", region = "Clarion-Clipperton Zone")
  )

  output <- assess_taxa_evidence_details(
    taxa,
    expected_environment = "marine",
    expected_habitat = "deep-sea benthic sediment",
    expected_habitat_type = "deep_sea_benthic",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )

  bird <- output[output$feature_id == "asv_bird", , drop = FALSE]
  reptile <- output[output$feature_id == "asv_reptile", , drop = FALSE]
  mammal <- output[output$feature_id == "asv_mammal", , drop = FALSE]
  expect_equal(bird$expected_habitat_status, "possible_unlikely")
  expect_equal(bird$ecological_status, "habitat_mismatch_or_allochthonous")
  expect_equal(bird$occurrence_interpretation, "possible_transient_or_allochthonous")
  expect_equal(bird$recommended_action, "flag_possible_exclusion")
  expect_equal(reptile$expected_habitat_status, "possible_unlikely")
  expect_equal(reptile$occurrence_interpretation, "possible_transient_or_allochthonous")
  expect_equal(reptile$recommended_action, "flag_possible_exclusion")
  expect_false(mammal$recommended_action == "flag_possible_exclusion")
  expect_false(mammal$expected_habitat_status == "possible_unlikely")
})

test_that("deep-sea benthic habitat type is inferred for nodule and seafloor wording", {
  taxa <- data.frame(
    feature_id = "asv_bird",
    phylum = "Chordata",
    class = "Aves",
    order = "Charadriiformes",
    family = "Laridae",
    genus = "Larus",
    species = "Larus dominicanus",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    "Larus dominicanus",
    "species",
    environment = "marine",
    region = "Clarion-Clipperton Zone"
  )

  output <- assess_taxa_evidence(
    taxa,
    expected_environment = "marine",
    expected_habitat = "abyssal polymetallic nodule seafloor",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )

  expect_equal(output$summary$recommended_action, "flag_possible_exclusion")
  expect_equal(output$summary$expected_habitat_status, "possible_unlikely")
  expect_equal(output$summary$rationale, "Marine seabird in deep-sea benthic sample; likely allochthonous.")
  expect_true(grepl("Aves or Reptilia", output$details$rationale, fixed = TRUE))
  expect_gt(nchar(output$details$rationale), nchar(output$summary$rationale))
})

test_that("plain deep-sea habitat infers deep-sea benthic bird rule", {
  taxa <- data.frame(
    feature_id = "asv_bird",
    phylum = "Chordata",
    class = "Aves",
    order = "Charadriiformes",
    family = "Laridae",
    genus = "Larus",
    species = "Larus dominicanus",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    "Larus dominicanus",
    "species",
    environment = "marine",
    region = "Cook Islands"
  )

  output <- assess_taxa_evidence(
    taxa,
    expected_environment = "marine",
    expected_habitat = "deep-sea",
    expected_region = "Cook Islands",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )

  expect_equal(output$summary$recommended_action, "flag_possible_exclusion")
  expect_equal(output$summary$expected_habitat_status, "possible_unlikely")
})

test_that("environment exclusion remains stronger than deep-sea bird/reptile rule", {
  taxa <- data.frame(
    feature_id = "asv_bird",
    phylum = "Chordata",
    class = "Aves",
    order = "Passeriformes",
    family = "Passeridae",
    genus = "Passer",
    species = "Passer domesticus",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    "Passer domesticus",
    "species",
    environment = "terrestrial",
    region = "Clarion-Clipperton Zone"
  )

  output <- assess_taxa_evidence(
    taxa,
    expected_environment = "marine",
    expected_habitat = "deep-sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct"
  )

  expect_equal(output$summary$expected_environment_status, "incompatible")
  expect_equal(output$summary$recommended_action, "exclude")
  expect_equal(output$summary$rationale, "Freshwater/terrestrial evidence conflicts with expected sample.")
})

test_that("deep-sea benthic pelagic taxa are retained as plausible transient signal", {
  tax <- data.frame(
    feature_id = "asv_skipjack",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Scombriformes",
    family = "Scombridae",
    genus = "Katsuwonus",
    species = "Katsuwonus pelamis",
    stringsAsFactors = FALSE
  )
  evidence <- assess_taxa_evidence_local_evidence(
    taxon_name = "Katsuwonus pelamis",
    taxon_rank = "species",
    environment = "marine",
    habitat = NA_character_,
    region = "Clarion-Clipperton Zone"
  )

  output <- assess_taxa_evidence_details(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep-sea benthic sediment",
    expected_habitat_type = "deep_sea_benthic",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    taxon_evidence_region_relation = "direct",
    ecology_evidence = assess_taxa_evidence_ecology_evidence()
  )

  expect_equal(output$expected_habitat_status, "possible_likely")
  expect_equal(output$recommended_action, "retain")
  expect_equal(output$occurrence_interpretation, "possible_transient_or_allochthonous")
})

test_that("list output works when source-specific evidence is absent", {
  output <- assess_taxa_evidence(
    assess_taxa_evidence_test_taxonomy()[1, , drop = FALSE],
    expected_environment = "marine",
    expected_habitat = "deep sea"
  )

  expect_true(all(c("summary", "details", "worms", "obis", "ecology") %in% names(output)))
  expect_s3_class(output$summary, "tbl_df")
  expect_s3_class(output$details, "data.frame")
  expect_true(inherits(output$obis, "data.frame"))
  expect_true(inherits(output$ecology, "data.frame"))
})
