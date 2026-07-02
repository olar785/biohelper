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
      recommended_action = c("retain", "flag_for_review"),
      rationale = c("Example rationale 1.", "Example rationale 2."),
      references = c("Example reference 1.", "Example reference 2."),
      stringsAsFactors = FALSE
    )
  )
}

valid_flag_taxa_structured_result <- function() {
  data.frame(
    feature_id = c("asv1", "asv2"),
    taxon_name = c("Salmo salar", "Daphnia"),
    taxon_rank = c("species", "genus"),
    expected_environment_status = c("compatible", "mixed_within_rank"),
    expected_habitat_status = c("transient_or_allochthonous_possible", "unknown"),
    expected_region_status = c("known_in_region", "known_near_region"),
    recommended_action = c("retain", "flag_for_review"),
    rationale = c("Example rationale 1.", "Example rationale 2."),
    references = c("Example reference 1.", "Example reference 2."),
    stringsAsFactors = FALSE
  )
}

flag_taxa_duplicate_taxonomy <- function() {
  data.frame(
    feature_id = c("asv1", "asv2", "asv3", "asv4", "asv5"),
    kingdom = rep("Animalia", 5),
    phylum = c("Chordata", "Chordata", "Arthropoda", "Porifera", NA_character_),
    class = c("Actinopteri", "Actinopteri", NA_character_, NA_character_, NA_character_),
    family = c("Salmonidae", "Salmonidae", NA_character_, NA_character_, NA_character_),
    genus = c("Salmo", "Salmo", "Daphnia", NA_character_, NA_character_),
    species = c("Salmo salar", "Salmo salar", NA_character_, NA_character_, NA_character_),
    raw_note = "raw taxonomy marker should not appear",
    stringsAsFactors = FALSE
  )
}

valid_flag_taxa_taxon_result <- function(expected_region_status = c("known_in_region", "unknown", "known_near_region")) {
  data.frame(
    taxon_name = c("Salmo salar", "Daphnia", "Porifera"),
    taxon_rank = c("species", "genus", "phylum"),
    expected_environment_status = c("compatible", "mixed_within_rank", "compatible"),
    expected_habitat_status = c("compatible", "unknown", "compatible"),
    expected_region_status = expected_region_status,
    recommended_action = c("retain", "flag_for_review", "retain"),
    rationale = c("Taxon-level rationale 1.", "Taxon-level rationale 2.", "Taxon-level rationale 3."),
    references = c("Taxon reference 1.", "Taxon reference 2.", "Taxon reference 3."),
    stringsAsFactors = FALSE
  )
}

valid_flag_taxa_taxon_result_with_tool_fields <- function(
  expected_region_status = c("known_in_region", "unknown", "known_near_region")
) {
  result <- valid_flag_taxa_taxon_result(expected_region_status = expected_region_status)
  result$evidence_basis <- c("local_evidence", "tool_evidence", "local_evidence")
  result$llm_tool_used <- c(FALSE, TRUE, FALSE)
  result$llm_tool_name <- c(NA_character_, "scite", NA_character_)
  result$llm_tool_query <- c(
    NA_character_,
    "Daphnia marine deep sea Clarion-Clipperton Zone",
    NA_character_
  )
  result$llm_tool_evidence_summary <- c(
    NA_character_,
    "Tool evidence found mixed marine and freshwater context for Daphnia.",
    NA_character_
  )
  result$llm_tool_references <- c(
    NA_character_,
    "Example Daphnia tool reference. doi:10.0000/example",
    NA_character_
  )
  result
}

add_query_ids <- function(result) {
  cbind(
    data.frame(
      query_id = sprintf("q%04d", seq_len(nrow(result))),
      stringsAsFactors = FALSE
    ),
    result,
    stringsAsFactors = FALSE
  )
}

flag_taxa_tool_evidence_result <- function(query_table) {
  data.frame(
    query_id = query_table$query_id,
    query_name = query_table$query_name,
    query_rank = query_table$query_rank,
    tool_used = TRUE,
    tool_name = "scite",
    tool_query = paste(query_table$query_name, "marine deep sea", sep = " "),
    tool_evidence_summary = paste("Explicit tool evidence for", query_table$query_name),
    tool_references = paste("Tool reference for", query_table$query_name),
    tool_found_explicit_evidence = TRUE,
    stringsAsFactors = FALSE
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

calanoida_like_taxonomy <- function(include_review_taxon = FALSE) {
  out <- data.frame(
    feature_id = "asv_calanoida",
    phylum = "Arthropoda",
    class = "Copepoda",
    order = "Calanoida",
    genus = NA_character_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_review_taxon)) {
    out <- rbind(
      out,
      data.frame(
        feature_id = "asv_daphnia",
        phylum = "Arthropoda",
        class = "Branchiopoda",
        order = NA_character_,
        genus = "Daphnia",
        stringsAsFactors = FALSE
      )
    )
  }
  out
}

calanoida_like_evidence <- function() {
  data.frame(
    taxon_name = c("Calanoida", "Calanoida"),
    taxon_rank = c("order", "order"),
    source = c("worms", "worms_ccz"),
    evidence_type = c("taxonomic_environment_database", "regional_deepsea_checklist"),
    evidence_summary = c(
      "WoRMS flags include marine and freshwater members.",
      "Calanoida is listed in the WoRMS Clarion-Clipperton Zone checklist."
    ),
    reference = c("WoRMS", "WoRMS CCZ checklist"),
    environment = c("marine; freshwater", "marine"),
    habitat = c(NA_character_, "deep sea"),
    region = c(NA_character_, "Clarion-Clipperton Zone"),
    reference_url = c(
      "https://www.marinespecies.org/aphia.php?p=taxdetails&id=1100",
      "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=1100"
    ),
    stringsAsFactors = FALSE
  )
}

prompt_allowed_value_block <- function(prompt, heading) {
  lines <- strsplit(prompt, "\n", fixed = TRUE)[[1]]
  start <- match(heading, lines)
  if (is.na(start)) {
    stop("Allowed-value heading not found: ", heading, call. = FALSE)
  }

  end <- start
  while (end < length(lines) && nzchar(lines[[end + 1]])) {
    end <- end + 1
  }

  paste(lines[start:end], collapse = "\n")
}

prompt_unique_query_table <- function(prompt) {
  lines <- strsplit(prompt, "\n", fixed = TRUE)[[1]]
  start <- match("Unique query taxa, tab-separated:", lines)
  if (is.na(start)) {
    stop("Unique query table marker not found.", call. = FALSE)
  }

  lines[(start + 1):length(lines)]
}

flag_taxa_prompt_for_test <- function(...) {
  args <- list(...)
  if (!("judgement_mode" %in% names(args))) {
    args$judgement_mode <- "llm_all"
  }
  result <- do.call(flag_taxa, args)
  if (is.character(result)) {
    return(result)
  }
  prompt_path <- attr(result, "prompt_path", exact = TRUE)
  if (is.null(prompt_path) || !file.exists(prompt_path)) {
    stop("flag_taxa prompt path was not available for test helper.", call. = FALSE)
  }
  paste(readLines(prompt_path, warn = FALSE), collapse = "\n")
}

flag_taxa_relevance_taxonomy <- function() {
  data.frame(
    feature_id = c("asv_species", "asv_family", "asv_unrelated"),
    kingdom = c("Animalia", "Animalia", "Animalia"),
    phylum = c("Chordata", "Porifera", "Mollusca"),
    class = c("Actinopteri", "Demospongiae", "Bivalvia"),
    family = c("Salmonidae", "Suberitidae", "Mytilidae"),
    genus = c("Salmo", NA_character_, "Mytilus"),
    species = c("Salmo salar", NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
}

flag_taxa_worms_evidence <- function() {
  data.frame(
    taxon_name = c("Salmo salar", "suberitidae", "Unrelated taxon"),
    taxon_rank = c("Species", "Family", "Species"),
    source = c("worms", "worms", "worms"),
    evidence_type = rep("taxonomic_environment_database", 3),
    evidence_summary = c(
      "Species evidence marker: Salmo salar appears marine-compatible in WoRMS.",
      "Family evidence marker: Suberitidae appears marine-compatible in WoRMS.",
      "Unrelated evidence marker that should not appear."
    ),
    reference = c(
      "WoRMS Editorial Board species",
      "WoRMS Editorial Board family",
      "WoRMS Editorial Board unrelated"
    ),
    accepted_name = c("Salmo salar", "Suberitidae", "Unrelated taxon"),
    accepted_rank = c("Species", "Family", "Species"),
    source_taxon_id = c("127186", "131628", "999999"),
    source_record_id = c("127186", "131628", "999999"),
    environment = c("marine", "marine", "terrestrial"),
    reference_url = c(
      "https://www.marinespecies.org/aphia.php?p=taxdetails&id=127186",
      "https://www.marinespecies.org/aphia.php?p=taxdetails&id=131628",
      "https://example.org/unrelated"
    ),
    stringsAsFactors = FALSE
  )
}

flag_taxa_ccz_evidence <- function() {
  data.frame(
    taxon_name = c("Salmo salar", "Unrelated CCZ taxon"),
    taxon_rank = c("Species", "Species"),
    source = c("worms_ccz", "worms_ccz"),
    evidence_type = rep("regional_deepsea_checklist", 2),
    evidence_summary = c(
      "CCZ regional marker: Salmo salar is listed in the CCZ checklist.",
      "Unrelated CCZ marker that should not appear."
    ),
    reference = c(
      "WoRMS CCZ checklist species reference",
      "WoRMS CCZ checklist unrelated reference"
    ),
    accepted_name = c("Salmo salar", "Unrelated CCZ taxon"),
    accepted_rank = c("Species", "Species"),
    source_taxon_id = c("127186", "888888"),
    source_record_id = c("127186", "888888"),
    environment = c("marine", "marine"),
    habitat = c("deep sea", "deep sea"),
    region = c("Clarion-Clipperton Zone", "Clarion-Clipperton Zone"),
    reference_url = c(
      "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=127186",
      "https://example.org/ccz/unrelated"
    ),
    kingdom = c("Animalia", "Animalia"),
    phylum = c("Chordata", "Mollusca"),
    class = c("Actinopteri", "Gastropoda"),
    order = c(NA_character_, NA_character_),
    family = c("Salmonidae", NA_character_),
    genus = c("Salmo", NA_character_),
    stringsAsFactors = FALSE
  )
}

flag_taxa_ccz_lineage_evidence <- function(n = 3) {
  data.frame(
    taxon_name = paste0("CCZ sponge ", seq_len(n)),
    taxon_rank = rep("Species", n),
    source = rep("worms_ccz", n),
    evidence_type = rep("regional_deepsea_checklist", n),
    evidence_summary = paste0("CCZ lineage marker ", seq_len(n), "."),
    reference = paste0("CCZ lineage reference ", seq_len(n)),
    accepted_name = paste0("Accepted CCZ sponge ", seq_len(n)),
    accepted_rank = rep("Species", n),
    source_taxon_id = as.character(9000 + seq_len(n)),
    source_record_id = as.character(9000 + seq_len(n)),
    environment = rep("marine", n),
    habitat = rep("deep sea", n),
    region = rep("Clarion-Clipperton Zone", n),
    reference_url = paste0("https://example.org/ccz/lineage/", seq_len(n)),
    kingdom = rep("Animalia", n),
    phylum = rep("Porifera", n),
    class = rep("Demospongiae", n),
    order = rep("Suberitida", n),
    family = rep("Suberitidae", n),
    genus = paste0("Cczgenus", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

flag_taxa_fetched_worms_evidence <- function() {
  data.frame(
    taxon_name = c("Salmo salar", "Daphnia"),
    taxon_rank = c("Species", "Genus"),
    source = c("worms", "worms"),
    evidence_type = rep("taxonomic_environment_database", 2),
    evidence_summary = c(
      "Fetched WoRMS evidence marker for Salmo salar.",
      "Fetched WoRMS evidence marker for Daphnia."
    ),
    reference = c("WoRMS fetched reference 1", "WoRMS fetched reference 2"),
    accepted_name = c("Salmo salar", "Daphnia"),
    source_taxon_id = c("127186", "12345"),
    source_record_id = c("127186", "12345"),
    environment = c("marine", "freshwater"),
    reference_url = c(
      "https://www.marinespecies.org/aphia.php?p=taxdetails&id=127186",
      "https://www.marinespecies.org/aphia.php?p=taxdetails&id=12345"
    ),
    stringsAsFactors = FALSE
  )
}

flag_taxa_load_ps_test_data_euk <- function() {
  testthat::skip_if_not_installed("phyloseq")
  data("ps_test_data_euk", package = "biohelper", envir = environment())
  ps_test_data_euk
}

test_that("prompt_only = TRUE returns a character prompt", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    family = "Salmonidae",
    genus = "Salmo",
    species = "Salmo salar"
  )

  prompt <- flag_taxa_prompt_for_test(
    tax,
    expected_environment = "marine",
    expected_habitat = "coastal water",
    expected_region = "North Atlantic"
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
  expect_true(grepl("Salmo salar", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE does not require ellmer or call the backend", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() {
      stop("ellmer should not be required.", call. = FALSE)
    },
    call_flag_taxa_ellmer = function(...) {
      stop("ellmer backend should not be called.", call. = FALSE)
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    prompt_only = TRUE
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
})

test_that("verbose reports unique taxa and dropped feature rows", {
  messages <- capture.output(
    prompt <- flag_taxa_prompt_for_test(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      verbose = TRUE
    ),
    type = "message"
  )

  expect_type(prompt, "character")
  expect_true(any(grepl("input contains 5 feature rows", messages)))
  expect_true(any(grepl("retained 4 feature rows with useful taxonomy; 1 will be returned as not assessed", messages)))
  expect_true(any(grepl("assessing 3 unique taxa from 4 retained feature rows", messages)))
  expect_true(any(grepl("Elapsed: [0-9.]+ sec", messages)))
})

test_that("verbose = FALSE suppresses diagnostic messages", {
  expect_silent(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      verbose = FALSE
    )
  )
})

test_that("prompt contains expected context values", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    genus = "Salmo",
    species = "Salmo salar"
  )

  prompt <- flag_taxa_prompt_for_test(
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

  prompt <- flag_taxa_prompt_for_test(
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
    "no_known_habitat_evidence",
    "transient_or_allochthonous_possible",
    "known_in_region",
    "known_near_region",
    "known_elsewhere_only",
    "restricted_elsewhere",
    "no_distribution_evidence",
    "not_assessed",
    "retain",
    "flag_for_review",
    "exclude"
  )

  for (allowed_value in allowed_values) {
    expect_true(grepl(allowed_value, prompt, fixed = TRUE))
  }
})

test_that("prompt excludes not_assessed from region status options when expected_region is provided", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone"
  )

  region_block <- prompt_allowed_value_block(prompt, "expected_region_status:")

  expect_true(grepl("- known_in_region", region_block, fixed = TRUE))
  expect_true(grepl("- no_distribution_evidence", region_block, fixed = TRUE))
  expect_true(grepl("- unknown", region_block, fixed = TRUE))
  expect_false(grepl("- not_assessed", region_block, fixed = TRUE))
  expect_true(grepl("Because expected_region is provided, do not use not_assessed.", prompt, fixed = TRUE))
  expect_true(grepl("If evidence is insufficient, use unknown or no_distribution_evidence as appropriate.", prompt, fixed = TRUE))
})

test_that("prompt includes not_assessed as the region status option when expected_region is NULL", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea"
  )

  region_block <- prompt_allowed_value_block(prompt, "expected_region_status:")

  expect_true(grepl("- not_assessed", region_block, fixed = TRUE))
  expect_false(grepl("- known_in_region", region_block, fixed = TRUE))
  expect_true(grepl("If expected_region is not specified, use not_assessed", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE asks for valid JSON only", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )

  json_fields <- c(
    "query_id",
    "taxon_name",
    "taxon_rank",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "recommended_action",
    "rationale",
    "references",
    "evidence_basis",
    "llm_tool_used",
    "llm_tool_name",
    "llm_tool_query",
    "llm_tool_evidence_summary",
    "llm_tool_references"
  )

  expect_true(grepl("Create a downloadable JSON file named `flag_taxa_llm_result.json`", prompt, fixed = TRUE))
  expect_true(grepl("output only the raw JSON array", prompt, fixed = TRUE))
  expect_true(grepl("valid JSON array of objects", prompt, fixed = TRUE))
  expect_true(grepl("Do not return markdown", prompt, fixed = TRUE))
  expect_true(grepl("plain text table", prompt, fixed = TRUE))
  expect_false(grepl("Return one table only", prompt, fixed = TRUE))

  for (json_field in json_fields) {
    expect_true(grepl(paste0("- ", json_field), prompt, fixed = TRUE))
  }
  expect_false(grepl("- feature_id", prompt, fixed = TRUE))
  expect_true(grepl("Do not include feature_id", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE uses unique query taxa instead of feature rows", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic"
  )
  query_table <- prompt_unique_query_table(prompt)
  query_header <- query_table[[1]]

  expect_equal(query_header, "query_id\tquery_name\tquery_rank\tlineage")
  expect_equal(sum(grepl("^q[0-9]+\\tSalmo salar\\tspecies\\t", query_table)), 1)
  expect_true(any(grepl("^q[0-9]+\\tDaphnia\\tgenus\\t", query_table)))
  expect_true(any(grepl("^q[0-9]+\\tPorifera\\tphylum\\t", query_table)))
  expect_false(any(grepl("^asv", query_table)))
  expect_false(grepl("feature_id", query_header, fixed = TRUE))
  expect_true(grepl("Phylum=Chordata", prompt, fixed = TRUE))
  expect_true(grepl("Class=Actinopteri", prompt, fixed = TRUE))
  expect_true(grepl("Order=NA", prompt, fixed = TRUE))
  expect_true(grepl("Family=Salmonidae", prompt, fixed = TRUE))
  expect_true(grepl("Genus=Salmo", prompt, fixed = TRUE))
  expect_true(grepl("Species=Salmo salar", prompt, fixed = TRUE))
  expect_true(grepl("Phylum=Porifera; Class=NA; Order=NA; Family=NA; Genus=NA; Species=NA", prompt, fixed = TRUE))
  expect_false(grepl("Kingdom=Animalia", prompt, fixed = TRUE))
  expect_false(grepl("raw taxonomy marker should not appear", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE can write one prompt to prompt_path", {
  prompt_file <- tempfile(fileext = ".txt")
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    prompt_path = prompt_file
  )

  expect_type(prompt, "character")
  expect_true(file.exists(prompt_file))
  expect_identical(paste(readLines(prompt_file, warn = FALSE), collapse = "\n"), prompt)
  expect_true(grepl("q0001\tSalmo salar", prompt, fixed = TRUE))
})

test_that("prompt writer stores all chunks in one txt file", {
  prompt_file <- tempfile(fileext = ".txt")
  prompts <- list(
    chunk_001 = "Required query_id values for this prompt: q0001\nq0001\tTaxon one",
    chunk_002 = "Required query_id values for this prompt: q0002\nq0002\tTaxon two"
  )

  written <- biohelper:::.write_flag_taxa_prompt_file(
    prompt_chunks = prompts,
    prompt_path = prompt_file
  )
  prompt <- paste(readLines(written, warn = FALSE), collapse = "\n")

  expect_equal(written, prompt_file)
  expect_true(file.exists(prompt_file))
  expect_true(grepl("===== CHUNK 1 / 2 =====", prompt, fixed = TRUE))
  expect_true(grepl("===== CHUNK 2 / 2 =====", prompt, fixed = TRUE))
  expect_true(grepl("q0001", prompt, fixed = TRUE))
  expect_true(grepl("q0002", prompt, fixed = TRUE))
  expect_equal(length(list.files(dirname(prompt_file), pattern = paste0("^", basename(prompt_file), "$"))), 1)
})

test_that("prompt_only reports prompt nchar when verbose", {
  messages <- capture.output(
    prompt <- flag_taxa_prompt_for_test(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = TRUE,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_type(prompt, "character")
  expect_true(any(grepl("prompt chunk 1 / 1 contains 3 taxa and [0-9]+ chars", messages)))
  expect_true(any(grepl("prompt includes 3 / 3 unique taxa selected for LLM/manual judgement", messages, fixed = TRUE)))
  expect_true(any(grepl("0 / 3 unique taxa were resolved deterministically and are not included in the prompt", messages, fixed = TRUE)))
})

test_that("unassessed feature rows are not included in prompt query table", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine"
  )
  query_table <- prompt_unique_query_table(prompt)

  expect_false(any(grepl("asv5", query_table, fixed = TRUE)))
  expect_false(any(grepl("\tNA\t", query_table, fixed = TRUE)))
})

test_that("allow_llm_tools controls prompt tool-use instructions", {
  prompt_without_tools <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    allow_llm_tools = FALSE
  )
  prompt_with_tools <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    judgement_mode = "llm_all",
    allow_llm_tools = TRUE,
    tool_requirement = "optional"
  )

  expect_true(grepl("Do not assume access to external tools or web search", prompt_without_tools, fixed = TRUE))
  expect_true(grepl("Use only the evidence provided in this prompt", prompt_without_tools, fixed = TRUE))
  expect_true(grepl("Registered tools may be available through the chat object", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Scite or another literature/search tool", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Use registered tools selectively", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("flagged and possible-exclusion boundary cases", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Do not use tools for taxa where local evidence clearly supports retain", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("local evidence is missing, weak, ambiguous, broad, contradictory, or indicates a possible resident-community concern", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Before assigning recommended_action = \"exclude\"", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Do not exclude a taxon based only on general model knowledge", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("If tool evidence does not clearly resolve the case, keep flag_for_review or flag_possible_exclusion as appropriate", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Use the full compact lineage to build precise search queries", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("Summarise any tool-derived evidence in rationale, references, and the tool-use fields", prompt_with_tools, fixed = TRUE))
  expect_true(grepl("final structured result", prompt_with_tools, fixed = TRUE))
})

test_that("tool_requirement default requires tool use for LLM-selected prompt rows", {
  expect_equal(eval(formals(flag_taxa)$tool_requirement), c("required_for_llm", "optional"))
  expect_equal(
    biohelper:::.validate_flag_taxa_tool_requirement("required_for_llm"),
    "required_for_llm"
  )
  expect_equal(
    biohelper:::.validate_flag_taxa_tool_requirement("optional"),
    "optional"
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    judgement_mode = "llm_all",
    allow_llm_tools = TRUE
  )

  expect_true(grepl("Tool use is mandatory for every query_id in this prompt", prompt, fixed = TRUE))
  expect_true(grepl("Use the registered literature/search tool available through the chat interface for every query_id in this prompt", prompt, fixed = TRUE))
  expect_true(grepl("llm_tool_references", prompt, fixed = TRUE))
  expect_true(grepl("Do not invent AphiaIDs, DOIs, URLs, article titles, author names, or citations", prompt, fixed = TRUE))
  expect_true(grepl("Use only references provided in the local evidence table or returned by tools", prompt, fixed = TRUE))
})

test_that("prompt_tools explicitly names requested manual prompt tools", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    judgement_mode = "llm_all",
    allow_llm_tools = TRUE,
    prompt_tools = "Scite"
  )

  expect_true(grepl("Use the Scite tool for every query_id in this prompt", prompt, fixed = TRUE))
  expect_true(grepl("For each query_id, report the tool query, tool evidence summary, and tool references", prompt, fixed = TRUE))
  expect_true(grepl("Do not invent references, AphiaIDs, DOIs, URLs, or citations", prompt, fixed = TRUE))
  expect_true(grepl("llm_tool_references", prompt, fixed = TRUE))
  expect_true(grepl("Do not use vague placeholders such as 'scientific_literature (tool evidence summary)'", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE does not run explicit tool pass by default", {
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(...) {
      stop("tool evidence pass should not be called.", call. = FALSE)
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    allow_llm_tools = TRUE,
    prompt_only = TRUE,
    chat = "fake_chat"
  )

  expect_type(prompt, "character")
  expect_false(grepl("Tool-derived evidence", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE does not run explicit tool pass", {
  calls <- 0
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(query_table, ...) {
      calls <<- calls + 1
      stop("tool evidence pass should not be called.", call. = FALSE)
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    judgement_mode = "llm_all",
    allow_llm_tools = TRUE,
    prompt_only = TRUE,
    chat = "fake_chat"
  )

  expect_equal(calls, 0)
  expect_type(prompt, "character")
  expect_false(grepl("Explicit tool evidence for Salmo salar", prompt, fixed = TRUE))
})

test_that("tool evidence pass respects max_tool_taxa and tool_batch_size", {
  batch_sizes <- integer()
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(query_table, ...) {
      batch_sizes <<- c(batch_sizes, nrow(query_table))
      flag_taxa_tool_evidence_result(query_table)
    },
    call_flag_taxa_ellmer = function(...) {
      add_query_ids(valid_flag_taxa_taxon_result_with_tool_fields(
        expected_region_status = rep("not_assessed", 3)
      ))
    }
  )

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    judgement_mode = "llm_all",
    allow_llm_tools = TRUE,
    max_tool_taxa = 2,
    tool_batch_size = 1,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_s3_class(output, "data.frame")
  expect_equal(batch_sizes, c(1L, 1L))
})

test_that("tool candidate selection follows judgement_mode", {
  tax_table <- extract_tax_table(flag_taxa_duplicate_taxonomy())
  feature_map <- biohelper:::.prepare_flag_taxa_feature_query_map(
    tax_table,
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
  )
  query_table <- biohelper:::.unique_flag_taxa_query_table(feature_map)
  preliminary <- data.frame(
    query_id = query_table$query_id,
    has_informative_local_evidence = c(TRUE, FALSE, TRUE),
    prelim_recommended_action = c("retain", "flag_for_review", "flag_possible_exclusion"),
    stringsAsFactors = FALSE
  )

  evidence_only <- biohelper:::.select_flag_taxa_tool_candidates(
    unique_query_table = query_table,
    preliminary_judgement = preliminary,
    judgement_mode = "evidence_only",
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = NULL,
    expected_region = NULL,
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
    max_evidence_rows_per_query = 5,
    tool_use_policy = "review_or_exclude"
  )
  expect_equal(nrow(evidence_only), 0)

  missing <- biohelper:::.select_flag_taxa_tool_candidates(
    unique_query_table = query_table,
    preliminary_judgement = preliminary,
    judgement_mode = "llm_missing_evidence",
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = NULL,
    expected_region = NULL,
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
    max_evidence_rows_per_query = 5,
    tool_use_policy = "review_or_exclude"
  )
  expect_equal(missing$query_id, query_table$query_id[[2]])

  flagged <- biohelper:::.select_flag_taxa_tool_candidates(
    unique_query_table = query_table,
    preliminary_judgement = preliminary,
    judgement_mode = "llm_flagged",
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
    max_evidence_rows_per_query = 5,
    tool_use_policy = "review_or_exclude"
  )
  expect_equal(flagged$query_id, query_table$query_id[2:3])

  possible_exclusion <- biohelper:::.select_flag_taxa_tool_candidates(
    unique_query_table = query_table,
    preliminary_judgement = preliminary,
    judgement_mode = "llm_possible_exclusion",
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
    max_evidence_rows_per_query = 5,
    tool_use_policy = "review_or_exclude"
  )
  expect_equal(possible_exclusion$query_id, query_table$query_id[[3]])

  review_alias <- biohelper:::.select_flag_taxa_tool_candidates(
    unique_query_table = query_table,
    preliminary_judgement = preliminary,
    judgement_mode = biohelper:::.validate_flag_taxa_judgement_mode(
      "llm_review",
      use_default = FALSE,
      prompt_only = FALSE,
      chat = "fake_chat",
      supplied_llm_result = NULL
    ),
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
    max_evidence_rows_per_query = 5,
    tool_use_policy = "review_or_exclude"
  )
  expect_equal(review_alias$query_id, flagged$query_id)

  all <- biohelper:::.select_flag_taxa_tool_candidates(
    unique_query_table = query_table,
    preliminary_judgement = preliminary,
    judgement_mode = "llm_all",
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = NULL,
    expected_region = NULL,
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
    max_evidence_rows_per_query = 5,
    tool_use_policy = "review_or_exclude"
  )
  expect_equal(all$query_id, query_table$query_id)
})

test_that("LLM query selection follows new judgement_mode values", {
  query_table <- data.frame(
    query_id = paste0("q000", 1:4),
    query_name = c("Retained", "Missing", "Review", "Possible exclusion"),
    query_rank = c("genus", "genus", "family", "class"),
    lineage = paste0("Phylum=Example; Class=Example; Order=NA; Family=NA; Genus=", 1:4, "; Species=NA"),
    stringsAsFactors = FALSE
  )
  preliminary <- data.frame(
    query_id = query_table$query_id,
    has_informative_local_evidence = c(TRUE, FALSE, TRUE, TRUE),
    prelim_recommended_action = c("retain", "flag_for_review", "flag_for_review", "flag_possible_exclusion"),
    stringsAsFactors = FALSE
  )

  expect_equal(
    biohelper:::.select_flag_taxa_llm_query_table(query_table, preliminary, "evidence_only")$query_id,
    character()
  )
  expect_equal(
    biohelper:::.select_flag_taxa_llm_query_table(query_table, preliminary, "llm_missing_evidence")$query_id,
    "q0002"
  )
  expect_equal(
    biohelper:::.select_flag_taxa_llm_query_table(query_table, preliminary, "llm_flagged")$query_id,
    c("q0002", "q0003", "q0004")
  )
  expect_equal(
    biohelper:::.select_flag_taxa_llm_query_table(query_table, preliminary, "llm_possible_exclusion")$query_id,
    "q0004"
  )
  expect_equal(
    biohelper:::.select_flag_taxa_llm_query_table(query_table, preliminary, "llm_all")$query_id,
    query_table$query_id
  )
  expect_equal(
    biohelper:::.validate_flag_taxa_judgement_mode(
      "llm_review",
      use_default = FALSE,
      prompt_only = FALSE,
      chat = "fake_chat",
      supplied_llm_result = NULL
    ),
    "llm_flagged"
  )
})

test_that("taxon_evidence = NULL keeps current prompt behaviour", {
  prompt_default <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )
  prompt_null <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    taxon_evidence = NULL
  )

  expect_identical(prompt_null, prompt_default)
})

test_that("valid taxon_evidence is accepted", {
  prompt <- flag_taxa_prompt_for_test(
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
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    taxon_evidence = valid_taxon_evidence()
  )

  expect_true(grepl("User-provided taxon evidence", prompt, fixed = TRUE))
  expect_true(grepl("Only evidence rows relevant to the input query taxa are included below.", prompt, fixed = TRUE))
  expect_true(grepl("Use this evidence before performing online searches.", prompt, fixed = TRUE))
  expect_true(grepl("Do not search online for taxa where the provided evidence is sufficient", prompt, fixed = TRUE))
  expect_true(grepl("Use online sources only when the provided evidence is incomplete, missing, or conflicting.", prompt, fixed = TRUE))
  expect_true(grepl("Compact relevant evidence table, tab-separated", prompt, fixed = TRUE))
  expect_true(grepl("Salmo salar", prompt, fixed = TRUE))
  expect_true(grepl("Example reference 1", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE does not include evidence section when taxon_evidence is NULL", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic"
  )

  expect_false(grepl("User-provided taxon evidence", prompt, fixed = TRUE))
})

test_that("evidence_sources = NULL does not fetch WoRMS evidence", {
  testthat::local_mocked_bindings(
    .flag_taxa_extract_taxa_for_evidence = function(...) {
      stop("extract_taxa_for_evidence should not be called.", call. = FALSE)
    },
    .flag_taxa_fetch_worms_evidence = function(...) {
      stop("fetch_worms_evidence should not be called.", call. = FALSE)
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    evidence_sources = NULL
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
  expect_false(grepl("Fetched WoRMS evidence marker", prompt, fixed = TRUE))
})

test_that("evidence_sources = worms fetches WoRMS evidence for prompt generation", {
  extract_calls <- 0
  fetch_calls <- 0

  testthat::local_mocked_bindings(
    .flag_taxa_extract_taxa_for_evidence = function(x) {
      extract_calls <<- extract_calls + 1
      expect_s3_class(x, "data.frame")
      data.frame(
        taxon_name = c("Salmo salar", "Daphnia"),
        taxon_rank = c("Species", "Genus"),
        stringsAsFactors = FALSE
      )
    },
    .flag_taxa_fetch_worms_evidence = function(taxa, by, cache_path, sleep, verbose, ...) {
      fetch_calls <<- fetch_calls + 1
      expect_equal(taxa$taxon_name, c("Salmo salar", "Daphnia"))
      expect_equal(by, "name")
      expect_null(cache_path)
      expect_equal(sleep, 0.2)
      expect_false(verbose)
      flag_taxa_fetched_worms_evidence()
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    evidence_sources = "worms"
  )

  expect_equal(extract_calls, 1)
  expect_equal(fetch_calls, 1)
  expect_true(grepl("Fetched WoRMS evidence marker for Salmo salar", prompt, fixed = TRUE))
  expect_true(grepl("WoRMS fetched reference 1", prompt, fixed = TRUE))
})

test_that("user and fetched WoRMS evidence are combined when both are supplied", {
  testthat::local_mocked_bindings(
    .flag_taxa_extract_taxa_for_evidence = function(x) {
      data.frame(
        taxon_name = c("Salmo salar", "Daphnia"),
        taxon_rank = c("Species", "Genus"),
        stringsAsFactors = FALSE
      )
    },
    .flag_taxa_fetch_worms_evidence = function(...) {
      flag_taxa_fetched_worms_evidence()
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = valid_taxon_evidence(),
    evidence_sources = "worms"
  )

  expect_true(grepl("Occurs in marine and freshwater phases.", prompt, fixed = TRUE))
  expect_true(grepl("Fetched WoRMS evidence marker for Salmo salar", prompt, fixed = TRUE))
  expect_true(grepl("Example reference 1", prompt, fixed = TRUE))
  expect_true(grepl("WoRMS fetched reference 1", prompt, fixed = TRUE))
})

test_that("unsupported evidence_sources errors clearly", {
  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      evidence_sources = "gbif"
    ),
    "unsupported value.*gbif.*worms"
  )
})

test_that("prompt_only = TRUE ignores mock_llm_result validation", {
  result <- valid_flag_taxa_result()
  result$recommended_action[1] <- "invalid"

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    prompt_only = TRUE,
    mock_llm_result = result
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
})

test_that("prompt_only = TRUE includes relevant exact species evidence", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "reef",
    expected_region = "Aotearoa New Zealand",
    taxon_evidence = flag_taxa_worms_evidence()
  )

  expect_true(grepl("Species evidence marker", prompt, fixed = TRUE))
  expect_true(grepl("Salmo salar", prompt, fixed = TRUE))
  expect_false(grepl("source_taxon_id", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE includes evidence for highest available higher-rank assignments", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = flag_taxa_worms_evidence()
  )

  expect_true(grepl("Family evidence marker", prompt, fixed = TRUE))
  expect_true(grepl("suberitidae", prompt, fixed = TRUE))
  expect_true(grepl("query_rank", prompt, fixed = TRUE))
  expect_true(grepl("family", prompt, fixed = TRUE))
})

test_that("prompt_only = TRUE excludes unrelated taxon evidence from relevant evidence section", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = flag_taxa_worms_evidence()
  )

  expect_false(grepl("Unrelated evidence marker", prompt, fixed = TRUE))
  expect_false(grepl("https://example.org/unrelated", prompt, fixed = TRUE))
})

test_that("accepted_name matches are included as relevant evidence", {
  evidence <- flag_taxa_worms_evidence()[1, , drop = FALSE]
  evidence$taxon_name <- "Synonym salmo"
  evidence$accepted_name <- "Salmo salar"
  evidence$evidence_summary <- "Accepted-name evidence marker."

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence
  )

  expect_true(grepl("Accepted-name evidence marker.", prompt, fixed = TRUE))
  expect_true(grepl("Synonym salmo", prompt, fixed = TRUE))
})

test_that("higher-rank lineage CCZ evidence is included for family queries", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = flag_taxa_ccz_lineage_evidence(3)
  )

  expect_true(grepl("CCZ lineage marker 1", prompt, fixed = TRUE))
  expect_true(grepl("Suberitidae", prompt, fixed = TRUE))
  expect_true(grepl("worms_ccz", prompt, fixed = TRUE))
  expect_true(grepl("regional_deepsea_checklist", prompt, fixed = TRUE))
  expect_true(grepl("Clarion-Clipperton Zone", prompt, fixed = TRUE))
  expect_true(grepl("deep sea", prompt, fixed = TRUE))
  expect_false(grepl("https://example.org/ccz/lineage/1", prompt, fixed = TRUE))
})

test_that("higher-rank lineage CCZ evidence is included for class and phylum queries", {
  tax <- data.frame(
    feature_id = c("asv_class", "asv_phylum"),
    kingdom = c("Animalia", "Animalia"),
    phylum = c("Porifera", "Porifera"),
    class = c("Demospongiae", NA_character_),
    order = c(NA_character_, NA_character_),
    family = c(NA_character_, NA_character_),
    genus = c(NA_character_, NA_character_),
    species = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )

  prompt <- flag_taxa_prompt_for_test(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = flag_taxa_ccz_lineage_evidence(3)
  )

  expect_true(grepl("query_rank\tquery_name", prompt, fixed = TRUE))
  expect_true(grepl("class\tDemospongiae", prompt, fixed = TRUE))
  expect_true(grepl("phylum\tPorifera", prompt, fixed = TRUE))
  expect_true(grepl("CCZ lineage marker 1", prompt, fixed = TRUE))
  expect_true(grepl("Demospongiae", prompt, fixed = TRUE))
  expect_true(grepl("Porifera", prompt, fixed = TRUE))
})

test_that("higher-rank evidence is not inferred for species-level query taxa", {
  tax <- data.frame(
    feature_id = "asv_species",
    kingdom = "Animalia",
    phylum = "Porifera",
    class = "Demospongiae",
    family = "Suberitidae",
    genus = "Aaptos",
    species = "Aaptos kanuux",
    stringsAsFactors = FALSE
  )
  evidence <- flag_taxa_worms_evidence()
  evidence <- evidence[evidence$taxon_name == "suberitidae", , drop = FALSE]

  prompt <- flag_taxa_prompt_for_test(
    tax,
    expected_environment = "marine",
    taxon_evidence = evidence
  )

  expect_true(grepl("Aaptos kanuux", prompt, fixed = TRUE))
  expect_false(grepl("Family evidence marker", prompt, fixed = TRUE))
})

test_that("WoRMS reference is preserved but reference_url is not pasted in the prompt", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = flag_taxa_worms_evidence()
  )

  expect_true(grepl("WoRMS Editorial Board species", prompt, fixed = TRUE))
  expect_false(grepl("https://www.marinespecies.org/aphia.php?p=taxdetails&id=127186", prompt, fixed = TRUE))
  expect_false(grepl("reference_url", prompt, fixed = TRUE))
})

test_that("long evidence references are compacted in the prompt", {
  evidence <- flag_taxa_worms_evidence()[1, , drop = FALSE]
  long_reference <- paste(rep("Very long WoRMS citation text", 20), collapse = " ")
  evidence$reference <- long_reference

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence
  )

  expect_false(grepl(long_reference, prompt, fixed = TRUE))
  expect_true(grepl("[truncated]", prompt, fixed = TRUE))
  expect_true(grepl("Very long WoRMS citation text", prompt, fixed = TRUE))
})

test_that("generic WoRMS and CCZ evidence for the same taxon are both retained", {
  evidence <- biohelper:::.combine_taxon_evidence_tables(
    flag_taxa_worms_evidence()[1, , drop = FALSE],
    flag_taxa_ccz_evidence()[1, , drop = FALSE]
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence
  )

  expect_true(grepl("Species evidence marker", prompt, fixed = TRUE))
  expect_true(grepl("CCZ regional marker", prompt, fixed = TRUE))
  expect_true(grepl("worms\t", prompt, fixed = TRUE))
  expect_true(grepl("worms_ccz", prompt, fixed = TRUE))
  expect_true(grepl("taxonomic_environment_database", prompt, fixed = TRUE))
  expect_true(grepl("regional_deepsea_checklist", prompt, fixed = TRUE))
  expect_true(grepl("deep sea", prompt, fixed = TRUE))
  expect_true(grepl("Clarion-Clipperton Zone", prompt, fixed = TRUE))
  expect_true(grepl("WoRMS CCZ checklist species reference", prompt, fixed = TRUE))
  expect_false(grepl(
    "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=127186",
    prompt,
    fixed = TRUE
  ))
})

test_that("generic WoRMS and CCZ evidence with different optional columns can both match one taxon", {
  tax <- data.frame(
    feature_id = "asv_plenaster",
    kingdom = "Animalia",
    phylum = "Porifera",
    class = "Demospongiae",
    order = "Tetractinellida",
    family = "Stellettidae",
    genus = "Plenaster",
    species = "Plenaster craigi",
    stringsAsFactors = FALSE
  )
  worms_evidence <- data.frame(
    taxon_name = "Plenaster craigi",
    taxon_rank = "Species",
    source = "worms",
    evidence_type = "taxonomic_environment_database",
    evidence_summary = "Generic WoRMS Plenaster marker.",
    reference = "WoRMS Plenaster reference",
    accepted_name = "Plenaster craigi",
    accepted_rank = "Species",
    source_taxon_id = "1590733",
    source_record_id = "1590733",
    environment = "marine",
    reference_url = "https://www.marinespecies.org/aphia.php?p=taxdetails&id=1590733",
    stringsAsFactors = FALSE
  )
  ccz_evidence <- data.frame(
    taxon_name = "Plenaster craigi",
    taxon_rank = "Species",
    source = "worms_ccz",
    evidence_type = "regional_deepsea_checklist",
    evidence_summary = "CCZ Plenaster marker.",
    reference = "WoRMS CCZ Plenaster reference",
    accepted_name = "Plenaster craigi",
    accepted_rank = "Species",
    source_taxon_id = "1590733",
    source_record_id = "1590733",
    environment = "marine",
    habitat = "deep sea",
    region = "Clarion-Clipperton Zone",
    reference_url = "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=1590733",
    kingdom = "Animalia",
    phylum = "Porifera",
    class = "Demospongiae",
    order = "Tetractinellida",
    family = "Stellettidae",
    genus = "Plenaster",
    stringsAsFactors = FALSE
  )
  combined_evidence <- dplyr::bind_rows(worms_evidence, ccz_evidence)

  prompt <- flag_taxa_prompt_for_test(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = combined_evidence,
    prompt_only = TRUE
  )

  expect_true(grepl("Generic WoRMS Plenaster marker.", prompt, fixed = TRUE))
  expect_true(grepl("CCZ Plenaster marker.", prompt, fixed = TRUE))
  expect_true(grepl("worms\t", prompt, fixed = TRUE))
  expect_true(grepl("worms_ccz", prompt, fixed = TRUE))
  expect_true(grepl("regional_deepsea_checklist", prompt, fixed = TRUE))
  expect_true(grepl("Clarion-Clipperton Zone", prompt, fixed = TRUE))
  expect_true(grepl("deep sea", prompt, fixed = TRUE))
  expect_false(grepl("https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=1590733", prompt, fixed = TRUE))
})

test_that("mocked chat succeeds when matched evidence has different optional columns", {
  tax <- data.frame(
    feature_id = "asv_plenaster",
    kingdom = "Animalia",
    phylum = "Porifera",
    family = "Stellettidae",
    genus = "Plenaster",
    species = "Plenaster craigi",
    stringsAsFactors = FALSE
  )
  combined_evidence <- dplyr::bind_rows(
    data.frame(
      taxon_name = "Plenaster craigi",
      taxon_rank = "Species",
      source = "worms",
      evidence_type = "taxonomic_environment_database",
      evidence_summary = "Generic WoRMS Plenaster marker.",
      reference = "WoRMS Plenaster reference",
      environment = "marine",
      stringsAsFactors = FALSE
    ),
    data.frame(
      taxon_name = "Plenaster craigi",
      taxon_rank = "Species",
      source = "worms_ccz",
      evidence_type = "regional_deepsea_checklist",
      evidence_summary = "CCZ Plenaster marker.",
      reference = "WoRMS CCZ Plenaster reference",
      environment = "marine",
      habitat = "deep sea",
      region = "Clarion-Clipperton Zone",
      reference_url = "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=1590733",
      phylum = "Porifera",
      family = "Stellettidae",
      genus = "Plenaster",
      stringsAsFactors = FALSE
    )
  )
  chat_result <- data.frame(
    feature_id = "asv_plenaster",
    taxon_name = "Plenaster craigi",
    taxon_rank = "species",
    expected_environment_status = "compatible",
    expected_habitat_status = "compatible",
    expected_region_status = "known_in_region",
    recommended_action = "retain",
    rationale = "Mocked LLM rationale.",
    references = "Mocked LLM reference.",
    stringsAsFactors = FALSE
  )
  calls <- 0
  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(
      prompt,
      chat,
      include_feature_id,
      expected_region,
      max_tries,
      retry_sleep
    ) {
      calls <<- calls + 1
      expect_true(grepl("Generic WoRMS Plenaster marker.", prompt, fixed = TRUE))
      expect_true(grepl("CCZ Plenaster marker.", prompt, fixed = TRUE))
      expect_equal(chat, "fake_chat")
      expect_false(include_feature_id)
      expect_equal(expected_region, "Clarion-Clipperton Zone")
      expect_equal(max_tries, 3)
      expect_equal(retry_sleep, 5)
      chat_result
    }
  )

  result <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = combined_evidence,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(calls, 1)
  expect_equal(result$feature_id, "asv_plenaster")
  expect_equal(result$taxon_name, "Plenaster craigi")
  expect_equal(result$expected_environment, "marine")
  expect_equal(result$expected_habitat, "deep sea")
  expect_equal(result$expected_region, "Clarion-Clipperton Zone")
  expect_true(grepl("Phylum=Porifera", result$lineage, fixed = TRUE))
  expect_true(grepl("worms", result$evidence_sources, fixed = TRUE))
  expect_true(grepl("worms_ccz", result$evidence_sources, fixed = TRUE))
  expect_true(grepl("environment=marine", result$evidence_summary, fixed = TRUE))
  expect_true(grepl("habitat=deep sea", result$evidence_summary, fixed = TRUE))
  expect_true(grepl("region=Clarion-Clipperton Zone", result$evidence_summary, fixed = TRUE))
})

test_that("unrelated CCZ evidence is excluded from the relevant evidence section", {
  evidence <- biohelper:::.combine_taxon_evidence_tables(
    flag_taxa_worms_evidence()[1, , drop = FALSE],
    flag_taxa_ccz_evidence()
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence
  )

  expect_true(grepl("CCZ regional marker", prompt, fixed = TRUE))
  expect_false(grepl("Unrelated CCZ marker", prompt, fixed = TRUE))
  expect_false(grepl("https://example.org/ccz/unrelated", prompt, fixed = TRUE))
})

test_that("unrelated rows from large CCZ evidence are not included in the prompt", {
  evidence <- biohelper:::.combine_taxon_evidence_tables(
    flag_taxa_ccz_lineage_evidence(2),
    data.frame(
      taxon_name = "Faraway CCZ species",
      taxon_rank = "Species",
      source = "worms_ccz",
      evidence_type = "regional_deepsea_checklist",
      evidence_summary = "Faraway marker that should not appear.",
      reference = "Faraway reference",
      accepted_name = "Faraway CCZ species",
      accepted_rank = "Species",
      environment = "marine",
      habitat = "deep sea",
      region = "Clarion-Clipperton Zone",
      reference_url = "https://example.org/ccz/faraway",
      phylum = "Cnidaria",
      class = "Anthozoa",
      family = "Farawayidae",
      stringsAsFactors = FALSE
    )
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence
  )

  expect_true(grepl("CCZ lineage marker 1", prompt, fixed = TRUE))
  expect_false(grepl("Faraway marker that should not appear.", prompt, fixed = TRUE))
  expect_false(grepl("https://example.org/ccz/faraway", prompt, fixed = TRUE))
})

test_that("many CCZ lineage matches are summarised rather than pasted row-by-row", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = flag_taxa_ccz_lineage_evidence(12)
  )

  expect_true(grepl("Summarised evidence:", prompt, fixed = TRUE))
  expect_true(grepl("matching evidence row", prompt, fixed = TRUE))
  expect_true(grepl("worms_ccz", prompt, fixed = TRUE))
  expect_true(grepl("regional_deepsea_checklist", prompt, fixed = TRUE))
  expect_true(grepl("Clarion-Clipperton Zone", prompt, fixed = TRUE))
  expect_true(grepl("deep sea", prompt, fixed = TRUE))
  expect_true(grepl("CCZ sponge 1", prompt, fixed = TRUE))
  expect_false(grepl("CCZ lineage marker 11", prompt, fixed = TRUE))
  expect_false(grepl("https://example.org/ccz/lineage/11", prompt, fixed = TRUE))
})

test_that("ccz_evidence_path combines local CCZ evidence with user evidence", {
  testthat::local_mocked_bindings(
    .flag_taxa_read_worms_ccz_evidence = function(path) {
      expect_equal(path, "CCZ_taxlist.csv")
      flag_taxa_ccz_evidence()[1, , drop = FALSE]
    }
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = flag_taxa_worms_evidence()[1, , drop = FALSE],
    ccz_evidence_path = "CCZ_taxlist.csv"
  )

  expect_true(grepl("Species evidence marker", prompt, fixed = TRUE))
  expect_true(grepl("CCZ regional marker", prompt, fixed = TRUE))
  expect_true(grepl("worms_ccz", prompt, fixed = TRUE))
})

test_that("WoRMS evidence is described as supporting evidence only", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = flag_taxa_worms_evidence()
  )

  expect_true(grepl("supporting evidence rather than an automatic decision", prompt, fixed = TRUE))
  expect_true(grepl("no-match, lookup failure, unknown, empty environment, or missing WoRMS evidence row is missing evidence, not evidence of incompatibility", prompt, fixed = TRUE))
})

test_that("flag_taxa does not call live WoRMS when evidence is supplied", {
  testthat::local_mocked_bindings(
    .require_worrms = function() stop("Live WoRMS should not be called.", call. = FALSE),
    .worms_records_names = function(...) stop("Live WoRMS should not be called.", call. = FALSE)
  )

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = flag_taxa_worms_evidence()
  )

  expect_type(prompt, "character")
  expect_length(prompt, 1)
})

test_that("flag_taxa includes relevant ps_test_data_euk evidence without live WoRMS calls", {
  ps <- flag_taxa_load_ps_test_data_euk()
  taxa <- extract_taxa_for_evidence(ps)
  taxa <- utils::head(taxa, 2)
  evidence <- data.frame(
    taxon_name = taxa$taxon_name,
    taxon_rank = taxa$taxon_rank,
    source = "worms",
    evidence_type = "taxonomic_environment_database",
    evidence_summary = paste("ps_test_data_euk evidence marker", taxa$taxon_name),
    reference = "WoRMS Editorial Board test",
    accepted_name = taxa$taxon_name,
    source_taxon_id = paste0("ps", seq_len(nrow(taxa))),
    source_record_id = paste0("ps", seq_len(nrow(taxa))),
    environment = "marine",
    reference_url = paste0("https://example.org/worms/", seq_len(nrow(taxa))),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    .require_worrms = function() stop("Live WoRMS should not be called.", call. = FALSE)
  )

  prompt <- flag_taxa_prompt_for_test(
    ps,
    expected_environment = "marine",
    taxon_evidence = evidence
  )

  expect_true(grepl("ps_test_data_euk evidence marker", prompt, fixed = TRUE))
  expect_true(grepl("WoRMS Editorial Board test", prompt, fixed = TRUE))
})

test_that("evidence_only works without chat or mock_llm_result", {
  tax <- data.frame(
    kingdom = "Animalia",
    phylum = "Chordata",
    species = "Salmo salar"
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    prompt_only = FALSE
  )

  expect_s3_class(output, "data.frame")
  expect_equal(output$recommended_action, "flag_for_review")
  expect_equal(output$evidence_basis, "conservative_reasoning_only")
})

test_that("max_prompt_chars errors before calling the LLM backend", {
  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(...) {
      stop("LLM backend should not be called.", call. = FALSE)
    }
  )

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      chat = "fake_chat",
      max_prompt_chars = 10
    ),
    "chunk 1.*1 taxa.*exceeds `max_prompt_chars.*Reduce `max_evidence_rows_per_query`.*max_tool_taxa.*larger-context model"
  )
})

test_that("local evidence summaries are truncated in prompt text only", {
  evidence <- flag_taxa_worms_evidence()[1, , drop = FALSE]
  long_summary <- paste(rep("very long local evidence sentence", 20), collapse = " ")
  evidence$evidence_summary <- long_summary

  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_relevance_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence
  )

  expect_true(grepl("[truncated]", prompt, fixed = TRUE))
  expect_false(grepl(long_summary, prompt, fixed = TRUE))
  expect_equal(evidence$evidence_summary, long_summary)
})

test_that("tool evidence summaries are truncated in internal prompt text only", {
  long_tool_summary <- paste(rep("explicit tool evidence sentence", 20), collapse = " ")
  query_table <- data.frame(
    query_id = "q0001",
    query_name = "Salmo salar",
    query_rank = "species",
    lineage = "Phylum=Chordata; Species=Salmo salar",
    stringsAsFactors = FALSE
  )
  tool_evidence <- flag_taxa_tool_evidence_result(query_table)
  tool_evidence$tool_evidence_summary <- long_tool_summary

  prompt <- biohelper:::.format_flag_taxa_tool_evidence_prompt_section(
    tool_evidence = tool_evidence,
    query_table = query_table,
    max_tool_evidence_summary_chars = 45
  )

  expect_true(grepl("[truncated]", prompt, fixed = TRUE))
  expect_false(grepl(long_tool_summary, prompt, fixed = TRUE))
})

test_that("lineage is truncated in internal prompt text only", {
  tax <- data.frame(
    phylum = "Chordata",
    class = "Actinopteri",
    order = "Salmoniformes",
    family = "Salmonidae",
    genus = "Salmo",
    species = "Salmo salar"
  )
  tax_table <- extract_tax_table(tax)
  feature_map <- biohelper:::.prepare_flag_taxa_feature_query_map(
    tax_table,
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
  )
  query_table <- biohelper:::.unique_flag_taxa_query_table(feature_map)
  preliminary <- biohelper:::.build_flag_taxa_preliminary_judgement(
    unique_query_table = query_table,
    taxon_evidence = NULL,
    expected_environment = "marine",
    expected_habitat = NULL,
    expected_region = NULL,
    tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
  )
  prompt <- biohelper:::build_flag_taxa_prompt(
    tax_table = tax_table,
    expected_environment = "marine",
    query_table = query_table,
    preliminary_judgement = preliminary,
    max_lineage_chars = 35
  )

  expect_true(grepl("Phylum=Chordata; Class=Actinopteri", prompt, fixed = TRUE))
  expect_true(grepl("[truncated]", prompt, fixed = TRUE))
})

test_that("auto_reduce_llm_chunk_size splits oversized chunks", {
  query_table <- data.frame(
    query_id = sprintf("q%04d", 1:4),
    query_name = paste0("taxon", 1:4),
    query_rank = "species",
    lineage = paste0("Phylum=P; Species=taxon", 1:4),
    stringsAsFactors = FALSE
  )
  chunks <- list(chunk_001 = query_table)
  reduced <- biohelper:::.auto_reduce_flag_taxa_query_chunks(
    query_chunks = chunks,
    build_prompt = function(chunk) paste(rep("x", nrow(chunk) * 1000), collapse = ""),
    max_prompt_chars = 2500
  )

  expect_equal(length(reduced), 2)
  expect_equal(unname(vapply(reduced, nrow, integer(1))), c(2L, 2L))
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

  expect_equal(nrow(output), nrow(tax))
  expect_equal(output$feature_id, tax$feature_id)
  expect_equal(output$taxon_name, c("Salmo salar", "Daphnia"))
  expect_equal(output$taxon_rank, c("species", "genus"))
  expect_equal(output$expected_environment, rep("marine", nrow(output)))
  expect_equal(output$expected_habitat, rep("estuary", nrow(output)))
  expect_equal(output$expected_region, rep("North Atlantic", nrow(output)))
  expect_true(all(grepl("Phylum=", output$lineage, fixed = TRUE)))
  expect_equal(output$recommended_action, result$recommended_action)
  expect_equal(output$rationale, result$rationale)
  expect_equal(output$references, result$references)
  expect_equal(output$evidence_sources, rep("none", nrow(output)))
  expect_equal(output$evidence_summary, rep("No matched local evidence.", nrow(output)))
})

test_that("valid compact mock_llm_result is returned", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_structured_result()

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_equal(output$feature_id, tax$feature_id)
  expect_equal(output$taxon_name, result$taxon_name)
  expect_equal(output$taxon_rank, result$taxon_rank)
  expect_equal(output$expected_environment, rep("marine", nrow(output)))
  expect_equal(output$expected_habitat, rep("estuary", nrow(output)))
  expect_equal(output$expected_region, rep("North Atlantic", nrow(output)))
  expect_equal(output$recommended_action, result$recommended_action)
  expect_equal(output$evidence_sources, rep("none", nrow(output)))
})

test_that("llm_result validates and joins external taxon-level results", {
  result <- add_query_ids(valid_flag_taxa_taxon_result_with_tool_fields())

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    llm_result = result
  )

  expect_equal(nrow(output), 5)
  expect_equal(output$feature_id, c("asv1", "asv2", "asv3", "asv4", "asv5"))
  expect_equal(output$taxon_name[1:2], rep("Salmo salar", 2))
  expect_equal(output$recommended_action[1:2], rep("retain", 2))
})

test_that("llm_result and mock_llm_result are mutually exclusive", {
  result <- valid_flag_taxa_taxon_result()

  expect_error(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      prompt_only = FALSE,
      llm_result = result,
      mock_llm_result = result
    ),
    "Provide only one of `llm_result` and `mock_llm_result`"
  )
})

test_that("taxon-level mock_llm_result is joined back to feature-level rows", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_equal(nrow(output), 5)
  expect_equal(output$feature_id, c("asv1", "asv2", "asv3", "asv4", "asv5"))
  expect_equal(output$taxon_name, c("Salmo salar", "Salmo salar", "Daphnia", "Porifera", NA_character_))
  expect_equal(output$recommended_action[1:2], c("retain", "retain"))
  expect_equal(output$recommended_action[[5]], "flag_for_review")
  expect_equal(output$expected_environment_status[[5]], "insufficient_taxonomic_resolution")
  expect_equal(
    output$rationale[[5]],
    "No useful taxonomic assignment was available, so this feature was not assessed by the LLM."
  )
  expect_equal(output$rationale[1:2], rep("Taxon-level rationale 1.", 2))
  expect_equal(output$evidence_basis[1:2], rep("local_evidence", 2))
  expect_equal(output$llm_tool_used[1:2], rep(FALSE, 2))
  expect_equal(output$evidence_basis[[3]], "tool_evidence")
  expect_true(output$llm_tool_used[[3]])
  expect_equal(output$llm_tool_name[[3]], "scite")
  expect_equal(
    output$llm_tool_query[[3]],
    "Daphnia marine deep sea Clarion-Clipperton Zone"
  )
  expect_equal(
    output$llm_tool_evidence_summary[[3]],
    "Tool evidence found mixed marine and freshwater context for Daphnia."
  )
  expect_equal(
    output$llm_tool_references[[3]],
    "Example Daphnia tool reference. doi:10.0000/example"
  )
  expect_equal(output$expected_environment, rep("marine", 5))
  expect_true(all(is.na(output$expected_habitat)))
  expect_equal(output$expected_region, rep("North Atlantic", 5))
  expect_true(grepl("Order=NA", output$lineage[[3]], fixed = TRUE))
  expect_equal(output$evidence_sources, rep("none", 5))
  expect_equal(output$evidence_summary[1:4], rep("No matched local evidence.", 4))
  expect_equal(output$evidence_summary[[5]], "No matched local evidence because no useful taxonomic assignment was available.")
})

test_that("features without useful taxonomy are returned as deterministic review rows", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  unassessed <- output[output$feature_id == "asv5", , drop = FALSE]
  expect_equal(nrow(output), nrow(flag_taxa_duplicate_taxonomy()))
  expect_equal(nrow(unassessed), 1)
  expect_true(is.na(unassessed$taxon_name))
  expect_true(is.na(unassessed$taxon_rank))
  expect_equal(unassessed$recommended_action, "flag_for_review")
  expect_equal(unassessed$expected_environment_status, "insufficient_taxonomic_resolution")
  expect_equal(unassessed$expected_habitat_status, "insufficient_taxonomic_resolution")
  expect_equal(unassessed$expected_region_status, "unknown")
  expect_true(grepl("not assessed by the LLM", unassessed$rationale, fixed = TRUE))
})

test_that("final output includes matched local evidence provenance", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()
  evidence <- biohelper:::.combine_taxon_evidence_tables(
    flag_taxa_worms_evidence()[1, , drop = FALSE],
    flag_taxa_ccz_evidence()[1, , drop = FALSE]
  )

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_true(all(c(
    "expected_environment",
    "expected_habitat",
    "expected_region",
    "lineage",
    "evidence_basis",
    "llm_tool_used",
    "llm_tool_name",
    "llm_tool_query",
    "llm_tool_evidence_summary",
    "llm_tool_references",
    "evidence_sources",
    "evidence_summary"
  ) %in% colnames(output)))
  expect_equal(output$evidence_sources[1:2], c("worms; worms_ccz", "worms; worms_ccz"))
  expect_true(grepl("worms: evidence_type=taxonomic_environment_database, environment=marine", output$evidence_summary[[1]], fixed = TRUE))
  expect_true(grepl("worms_ccz: evidence_type=regional_deepsea_checklist, environment=marine, habitat=deep sea, region=Clarion-Clipperton Zone", output$evidence_summary[[1]], fixed = TRUE))
  expect_equal(output$evidence_sources[3:4], c("none", "none"))
  expect_equal(output$evidence_summary[3:4], rep("No matched local evidence.", 2))
})

test_that("old taxon-level mock_llm_result gets conservative tool defaults", {
  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    mock_llm_result = valid_flag_taxa_taxon_result()
  )

  expect_equal(output$evidence_basis, rep("conservative_reasoning_only", nrow(output)))
  expect_equal(output$llm_tool_used, rep(FALSE, nrow(output)))
  expect_true(all(is.na(output$llm_tool_name)))
  expect_true(all(is.na(output$llm_tool_query)))
  expect_true(all(is.na(output$llm_tool_evidence_summary)))
  expect_true(all(is.na(output$llm_tool_references)))
})

test_that("taxon-level mock_llm_result with not_assessed passes when expected_region is NULL", {
  result <- valid_flag_taxa_taxon_result(expected_region_status = "not_assessed")

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_equal(unique(output$expected_region_status), "not_assessed")
  expect_true(all(is.na(output$expected_region)))
})

test_that("invalid taxon-level mock_llm_result fails validation", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()
  result$recommended_action[1] <- "maybe"

  expect_error(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "recommended_action.*invalid values.*maybe"
  )
})

test_that("invalid evidence_basis fails validation", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()
  result$evidence_basis[1] <- "general_knowledge"

  expect_error(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = result
    ),
    "evidence_basis.*invalid values.*general_knowledge"
  )
})

test_that("llm_tool_used = TRUE requires a tool name and evidence summary", {
  missing_name <- valid_flag_taxa_taxon_result_with_tool_fields()
  missing_name$llm_tool_name[2] <- NA_character_

  expect_error(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = missing_name
    ),
    "llm_tool_name.*must not be missing or empty"
  )

  missing_summary <- valid_flag_taxa_taxon_result_with_tool_fields()
  missing_summary$llm_tool_evidence_summary[2] <- ""

  expect_error(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = missing_summary
    ),
    "llm_tool_evidence_summary.*must not be missing or empty"
  )
})

test_that("missing llm_tool_references are repaired by default", {
  missing_refs <- valid_flag_taxa_taxon_result_with_tool_fields()
  missing_refs$llm_tool_references[2] <- ""

  expect_message(
    output <- flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      verbose = TRUE,
      mock_llm_result = missing_refs
    ),
    "filled missing llm_tool_references for 1 tool-used row"
  )

  expect_equal(
    output$llm_tool_references[[3]],
    "Tool used but no verifiable references were returned."
  )
  expect_true(all(is.na(output$llm_tool_references[c(1, 2, 4, 5)])))

  real_refs <- valid_flag_taxa_taxon_result_with_tool_fields()
  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    mock_llm_result = real_refs
  )
  expect_equal(
    output$llm_tool_references[[3]],
    "Example Daphnia tool reference. doi:10.0000/example"
  )
})

test_that("required tool attempts are repaired conservatively by default", {
  result <- valid_flag_taxa_taxon_result()
  result$recommended_action[2] <- "exclude"

  expect_message(
    output <- flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      allow_llm_tools = TRUE,
      prompt_only = FALSE,
      verbose = TRUE,
      mock_llm_result = result
    ),
    "filled required tool-evidence metadata for 3 LLM-selected rows"
  )

  expect_true(all(output$llm_tool_used[1:4]))
  expect_equal(unique(output$llm_tool_name[1:4]), "registered_tool")
  expect_true(all(output$llm_tool_references[1:4] == "Tool evidence was required but no verifiable references were returned."))
  expect_true(all(grepl("Tool evidence was required but no usable tool evidence was returned", output$llm_tool_evidence_summary[1:4], fixed = TRUE)))
  expect_equal(output$recommended_action[[3]], "flag_for_review")
  expect_true(grepl("biohelper repair: tool evidence was required", output$rationale[[3]], fixed = TRUE))
})

test_that("strict required-tool validation can still error", {
  query_table <- biohelper:::.unique_flag_taxa_query_table(
    biohelper:::.prepare_flag_taxa_feature_query_map(
      extract_tax_table(flag_taxa_duplicate_taxonomy()),
      tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
    )
  )

  expect_error(
    biohelper:::.repair_flag_taxa_required_tool_attempts(
      result = add_query_ids(valid_flag_taxa_taxon_result()),
      tool_requirement = "required_for_llm",
      allow_llm_tools = TRUE,
      repair_invalid_llm = FALSE,
      query_taxonomy = query_table
    ),
    "requires every LLM-selected row to report a tool attempt"
  )
})

test_that("strict tool reference validation still errors without repair", {
  query_table <- biohelper:::.unique_flag_taxa_query_table(
    biohelper:::.prepare_flag_taxa_feature_query_map(
      extract_tax_table(flag_taxa_duplicate_taxonomy()),
      tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
    )
  )
  missing_refs <- add_query_ids(valid_flag_taxa_taxon_result_with_tool_fields())
  missing_refs$llm_tool_references[2] <- ""
  unrepaired <- biohelper:::.repair_flag_taxa_invalid_llm_output(
    result = missing_refs,
    repair_invalid_llm = FALSE,
    expected_region = "North Atlantic"
  )

  expect_error(
    biohelper:::validate_flag_taxa_taxon_output(
      result = unrepaired,
      query_taxonomy = query_table,
      expected_region = "North Atlantic",
      allow_llm_tools = TRUE
    ),
    "llm_tool_references.*must not be missing or empty"
  )
})

test_that("vague tool references still error and explicit no-reference text passes", {
  vague_refs <- valid_flag_taxa_taxon_result_with_tool_fields()
  vague_refs$llm_tool_references[2] <- "scientific_literature (tool evidence summary)"

  expect_error(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = vague_refs
    ),
    "Vague placeholders"
  )

  no_refs <- valid_flag_taxa_taxon_result_with_tool_fields()
  no_refs$llm_tool_references[2] <- "Tool used but no verifiable references were returned."
  expect_s3_class(
    flag_taxa(
      flag_taxa_duplicate_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      mock_llm_result = no_refs
    ),
    "data.frame"
  )
})

test_that("conservative reasoning cannot support exclude", {
  result <- add_query_ids(valid_flag_taxa_taxon_result_with_tool_fields())
  result$recommended_action[1] <- "exclude"
  result$evidence_basis[1] <- "conservative_reasoning_only"
  result$llm_tool_used[1] <- FALSE
  query_table <- biohelper:::.unique_flag_taxa_query_table(
    biohelper:::.prepare_flag_taxa_feature_query_map(
      extract_tax_table(flag_taxa_duplicate_taxonomy()),
      tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
    )
  )

  expect_error(
    biohelper:::validate_flag_taxa_taxon_output(
      result = result,
      query_taxonomy = query_table,
      expected_region = "North Atlantic",
      allow_llm_tools = TRUE
    ),
    "exclude.*conservative_reasoning_only"
  )
})

test_that("review with conservative reasoning passes and retain with local evidence passes", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()
  result$recommended_action[1] <- "flag_for_review"
  result$evidence_basis[1] <- "conservative_reasoning_only"
  result$llm_tool_used[1] <- FALSE
  result$llm_tool_name[1] <- NA_character_
  result$llm_tool_query[1] <- NA_character_
  result$llm_tool_evidence_summary[1] <- NA_character_
  result$recommended_action[3] <- "retain"
  result$evidence_basis[3] <- "local_evidence"
  result$llm_tool_used[3] <- FALSE

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    allow_llm_tools = TRUE,
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_equal(output$recommended_action[[1]], "flag_for_review")
  expect_equal(output$evidence_basis[[1]], "conservative_reasoning_only")
  expect_equal(output$recommended_action[[4]], "retain")
  expect_equal(output$evidence_basis[[4]], "local_evidence")
})

test_that("prompt_only = FALSE with a mocked chat backend returns a validated result", {
  calls <- 0
  result <- valid_flag_taxa_structured_result()

  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(
      prompt,
      chat,
      include_feature_id,
      expected_region,
      max_tries,
      retry_sleep
    ) {
      calls <<- calls + 1
      expect_type(prompt, "character")
      expect_equal(chat, "fake_chat")
      expect_false(include_feature_id)
      expect_equal(expected_region, "North Atlantic")
      expect_equal(max_tries, 3)
      expect_equal(retry_sleep, 5)
      result
    }
  )

  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(calls, 1)
  expect_equal(output$feature_id, result$feature_id)
  expect_equal(output$taxon_name, result$taxon_name)
  expect_equal(output$taxon_rank, result$taxon_rank)
  expect_equal(output$expected_environment, rep("marine", nrow(output)))
  expect_equal(output$expected_habitat, rep("estuary", nrow(output)))
  expect_equal(output$expected_region, rep("North Atlantic", nrow(output)))
  expect_equal(output$recommended_action, result$recommended_action)
})

test_that("fake chat returning taxon-level output is joined back to features", {
  calls <- 0
  result <- valid_flag_taxa_taxon_result_with_tool_fields()

  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(
      prompt,
      chat,
      include_feature_id,
      expected_region,
      max_tries,
      retry_sleep
    ) {
      calls <<- calls + 1
      expect_type(prompt, "character")
      expect_equal(chat, "fake_chat")
      expect_false(include_feature_id)
      expect_equal(expected_region, "North Atlantic")
      result
    }
  )

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(calls, 1)
  expect_equal(nrow(output), 5)
  expect_equal(output$feature_id, c("asv1", "asv2", "asv3", "asv4", "asv5"))
  expect_equal(output$taxon_name[1:2], c("Salmo salar", "Salmo salar"))
  expect_equal(output$rationale[1:2], rep("Taxon-level rationale 1.", 2))
  expect_equal(output$recommended_action[1:2], rep("retain", 2))
  expect_equal(output$llm_tool_name[[3]], "scite")
  expect_true(is.na(output$taxon_name[[5]]))
})

test_that("llm_review output categorical columns are string labels, not enum codes", {
  calls <- 0

  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(
      prompt,
      chat,
      include_feature_id,
      expected_region,
      max_tries,
      retry_sleep
    ) {
      calls <<- calls + 1
      expect_equal(chat, "fake_chat")
      data.frame(
        query_id = sprintf("q%04d", 1:3),
        taxon_name = c("Salmo salar", "Daphnia", "Porifera"),
        taxon_rank = c("species", "genus", "phylum"),
        expected_environment_status = c(1L, 3L, 1L),
        expected_habitat_status = c(1L, 8L, 1L),
        expected_region_status = c(1L, 6L, 2L),
        recommended_action = c(1L, 2L, 1L),
        rationale = c("Numeric enum 1.", "Numeric enum 2.", "Numeric enum 3."),
        references = c("Reference 1.", "Reference 2.", "Reference 3."),
        evidence_basis = c(4L, 4L, 4L),
        llm_tool_used = FALSE,
        llm_tool_name = NA_character_,
        llm_tool_query = NA_character_,
        llm_tool_evidence_summary = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  )

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    judgement_mode = "llm_review",
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(calls, 1)
  expect_equal(nrow(output), 5)
  categorical_columns <- c(
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "recommended_action",
    "evidence_basis"
  )
  expect_true(all(vapply(output[categorical_columns], is.character, logical(1))))
  expect_false(any(vapply(
    output[categorical_columns],
    function(column) any(grepl("^[0-9]+$", column)),
    logical(1)
  )))
  expect_equal(output$expected_environment_status[[1]], "compatible")
  expect_equal(output$expected_environment_status[[3]], "mixed_within_rank")
  expect_equal(output$expected_habitat_status[[3]], "unknown")
  expect_equal(output$recommended_action[[3]], "flag_for_review")
})

test_that("tool evidence pass is included in final prompt and propagated to result", {
  tool_calls <- 0
  llm_calls <- 0
  result <- valid_flag_taxa_taxon_result()

  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(query_table, ...) {
      tool_calls <<- tool_calls + 1
      flag_taxa_tool_evidence_result(query_table)
    },
    call_flag_taxa_ellmer = function(
      prompt,
      chat,
      include_feature_id,
      expected_region,
      max_tries,
      retry_sleep
    ) {
      llm_calls <<- llm_calls + 1
      expect_true(grepl("Tool-derived evidence", prompt, fixed = TRUE))
      expect_true(grepl("Explicit tool evidence for Salmo salar", prompt, fixed = TRUE))
      result
    }
  )

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_region = "North Atlantic",
    allow_llm_tools = TRUE,
    tool_requirement = "optional",
    judgement_mode = "llm_all",
    max_tool_taxa = 1,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(tool_calls, 1)
  expect_equal(llm_calls, 1)
  expect_true(output$llm_tool_used[[1]])
  expect_equal(output$llm_tool_name[[1]], "scite")
  expect_equal(output$evidence_basis[[1]], "tool_evidence")
  expect_true(grepl("Explicit tool evidence for Salmo salar", output$llm_tool_evidence_summary[[1]], fixed = TRUE))
  expect_false(output$llm_tool_used[[3]])
})

test_that("required tool pass failure is recorded as conservative tool metadata", {
  evidence <- valid_taxon_evidence()[1, , drop = FALSE]
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(...) {
      stop("simulated tool outage", call. = FALSE)
    },
    call_flag_taxa_ellmer = function(prompt, ...) {
      expect_true(grepl("tool_failed:", prompt, fixed = TRUE))
      data.frame(
        query_id = "q0002",
        taxon_name = "Daphnia",
        taxon_rank = "genus",
        expected_environment_status = "unknown",
        expected_habitat_status = "no_known_habitat_evidence",
        expected_region_status = "not_assessed",
        ecological_status = "unknown",
        occurrence_interpretation = "uncertain",
        recommended_action = "flag_for_review",
        rationale = "LLM remained conservative after tool failure.",
        references = "No verifiable tool references.",
        evidence_basis = "conservative_reasoning_only",
        llm_tool_used = FALSE,
        llm_tool_name = NA_character_,
        llm_tool_query = NA_character_,
        llm_tool_evidence_summary = NA_character_,
        llm_tool_references = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  )

  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence,
    judgement_mode = "llm_missing_evidence",
    allow_llm_tools = TRUE,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_true(output$llm_tool_used[[2]])
  expect_equal(output$llm_tool_name[[2]], "registered_tool")
  expect_true(grepl("tool_failed:", output$llm_tool_evidence_summary[[2]], fixed = TRUE))
  expect_equal(
    output$llm_tool_references[[2]],
    "Tool evidence was required but no verifiable references were returned."
  )
  expect_equal(output$recommended_action[[2]], "flag_for_review")
})

test_that("prompt says no explicit evidence is uncertainty, not exclusion", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    allow_llm_tools = TRUE,
    prompt_only = TRUE
  )

  expect_true(grepl("Failure to find explicit evidence", prompt, fixed = TRUE))
  expect_true(grepl("No evidence found' means uncertainty, not exclusion", prompt, fixed = TRUE))
  expect_true(grepl("no_known_habitat_evidence", prompt, fixed = TRUE))
  expect_true(grepl("no Scite evidence", prompt, fixed = TRUE))
})

test_that("no_known_habitat_evidence is allowed for review outputs", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()
  result$expected_habitat_status[[2]] <- "no_known_habitat_evidence"
  result$recommended_action[[2]] <- "flag_for_review"
  result$llm_tool_evidence_summary[[2]] <- "no_evidence_found: Tool search found no explicit evidence for this taxon in the expected habitat/region; this is missing evidence, not incompatibility."
  result$references[[2]] <- "Tool search returned no explicit evidence."

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    prompt_only = FALSE,
    llm_result = result
  )

  expect_equal(output$expected_habitat_status[[3]], "no_known_habitat_evidence")
  expect_equal(output$recommended_action[[3]], "flag_for_review")
})

test_that("evidence_only keeps broad marine-supported habitat gaps as review", {
  tax <- data.frame(
    feature_id = "asv_arcellinida",
    phylum = "Amoebozoa",
    class = "Tubulinea",
    order = "Arcellinida",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Arcellinida",
    taxon_rank = "order",
    source = "worms",
    evidence_type = "taxonomic_environment_database",
    evidence_summary = "WoRMS environment flags include marine; brackish; freshwater; terrestrial.",
    reference = "WoRMS",
    environment = "marine; brackish; freshwater; terrestrial",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$expected_environment_status, "mixed_within_rank")
  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_equal(output$recommended_action, "flag_for_review")
  expect_false(output$recommended_action == "exclude")
})

test_that("Calanoida-like mixed environment with compatible habitat and region is retained", {
  output <- flag_taxa(
    calanoida_like_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = calanoida_like_evidence(),
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$expected_environment_status, "mixed_within_rank")
  expect_equal(output$expected_habitat_status, "compatible")
  expect_equal(output$expected_region_status, "known_in_region")
  expect_equal(output$ecological_status, "compatible")
  expect_equal(output$occurrence_interpretation, "expected_resident")
  expect_equal(output$recommended_action, "retain")
  expect_equal(output$worms_environment, "marine; freshwater")
})

test_that("Arcellinida-like questionable broad marine evidence is possible exclusion, not hard exclude", {
  tax <- data.frame(
    feature_id = "asv_arcellinida",
    phylum = "Amoebozoa",
    class = "Tubulinea",
    order = "Arcellinida",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Arcellinida",
    taxon_rank = "order",
    source = "worms",
    evidence_type = "taxonomic_environment_database",
    evidence_summary = "Most true Arcellinida are freshwater or terrestrial; marine records may be questionable.",
    reference = "Example Arcellinida ecology",
    environment = "marine; freshwater; terrestrial",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$expected_environment_status, "mixed_within_rank")
  expect_equal(output$ecological_status, "unlikely_resident")
  expect_equal(output$occurrence_interpretation, "possible_contaminant_or_misassignment")
  expect_equal(output$recommended_action, "flag_possible_exclusion")
  expect_true(grepl("predominantly non-marine", output$rationale, fixed = TRUE))
  expect_true(grepl("positive ecological concern", output$rationale, fixed = TRUE))
  expect_false(output$recommended_action == "exclude")
})

test_that("Oomycota-like host-associated possible exclusion explains resident-habitat concern", {
  tax <- data.frame(
    feature_id = "asv_oomycota",
    phylum = "Oomycota",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Oomycota",
    taxon_rank = "phylum",
    source = "literature",
    evidence_type = "ecology",
    evidence_summary = "Oomycota includes host-associated parasitic and pathogenic taxa with marine members.",
    reference = "Oomycota ecology",
    environment = "marine",
    habitat = "host-associated; parasitic",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$recommended_action, "flag_possible_exclusion")
  expect_equal(output$ecological_status, "unlikely_resident")
  expect_true(grepl("host-associated", output$rationale, fixed = TRUE))
  expect_true(grepl("free-living resident", output$rationale, fixed = TRUE))
})

test_that("Discosea-like broad non-marine ecology possible exclusion explains ecological concern", {
  tax <- data.frame(
    feature_id = "asv_discosea",
    phylum = "Amoebozoa",
    class = "Discosea",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Discosea",
    taxon_rank = "class",
    source = "literature",
    evidence_type = "ecology",
    evidence_summary = "Discosea are primarily freshwater, soil, and terrestrial amoebae; marine records may be questionable.",
    reference = "Discosea ecology",
    environment = "marine; freshwater; terrestrial",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$recommended_action, "flag_possible_exclusion")
  expect_equal(output$ecological_status, "unlikely_resident")
  expect_true(grepl("questionable marine records", output$rationale, fixed = TRUE))
  expect_true(grepl("resident marine/deep-sea interpretation", output$rationale, fixed = TRUE))
})

test_that("photosynthetic Florideophyceae in aphotic deep sea is resident-incompatible but allochthonous possible", {
  tax <- data.frame(
    feature_id = "asv_florideophyceae",
    phylum = "Rhodophyta",
    class = "Florideophyceae",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Florideophyceae",
    taxon_rank = "class",
    source = "literature",
    evidence_type = "ecology",
    evidence_summary = "Florideophyceae is a photosynthetic red algal class with marine members.",
    reference = "Example red algal ecology",
    environment = "marine",
    habitat = "photic coastal water",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "aphotic deep-sea benthic sediment",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$expected_environment_status, "compatible")
  expect_equal(output$expected_habitat_status, "incompatible")
  expect_equal(output$ecological_status, "incompatible_resident")
  expect_equal(output$occurrence_interpretation, "possible_transient_or_allochthonous")
  expect_true(output$recommended_action %in% c("flag_possible_exclusion", "exclude"))
  expect_true(grepl("photosynthetic", output$rationale, fixed = TRUE))
  expect_true(grepl("aphotic/deep-sea/benthic", output$rationale, fixed = TRUE))
  expect_true(grepl("allochthonous DNA remains possible", output$rationale, fixed = TRUE))
})

test_that("prompt_only returns evidence-first result and writes one txt prompt by default", {
  old_wd <- getwd()
  temp_wd <- tempfile("flag_taxa_prompt_wd_")
  dir.create(temp_wd)
  setwd(temp_wd)
  on.exit(setwd(old_wd), add = TRUE)

  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    judgement_mode = "llm_all",
    prompt_only = TRUE
  )

  prompt_path <- attr(output, "prompt_path", exact = TRUE)
  expect_s3_class(output, "data.frame")
  expect_true(file.exists(prompt_path))
  expect_equal(basename(prompt_path), "flag_taxa_prompt.txt")
  prompt <- paste(readLines(prompt_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("Create a downloadable JSON file named `flag_taxa_llm_result.json`", prompt, fixed = TRUE))
  expect_true(grepl("Required query_id values for this prompt", prompt, fixed = TRUE))
  expect_true(all(c("ecological_status", "occurrence_interpretation", "recommended_action") %in% colnames(output)))
  expect_true("llm_prompt_selected" %in% colnames(output))
})

test_that("prompt_only output marks prompt-selected taxa and written prompt includes all selected IDs", {
  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = valid_taxon_evidence()[1, , drop = FALSE],
    judgement_mode = "llm_missing_evidence",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  prompt <- paste(readLines(attr(output, "prompt_path", exact = TRUE), warn = FALSE), collapse = "\n")
  selected_ids <- unique(as.character(output$query_id[output$llm_prompt_selected]))
  selected_ids <- selected_ids[!is.na(selected_ids)]
  bypass_ids <- unique(as.character(output$query_id[!output$llm_prompt_selected]))
  bypass_ids <- bypass_ids[!is.na(bypass_ids)]

  expect_true(length(selected_ids) > 0)
  expect_equal(output$llm_prompt_selected, output$llm_selected)
  for (query_id in selected_ids) {
    expect_true(grepl(query_id, prompt, fixed = TRUE))
  }
  for (query_id in bypass_ids) {
    expect_false(grepl(paste0(query_id, "\t"), prompt, fixed = TRUE))
  }
})

test_that("prompt_only selection follows judgement_mode", {
  tax <- data.frame(
    feature_id = c("asv_calanoida", "asv_oomycota", "asv_daphnia"),
    phylum = c("Arthropoda", "Oomycota", "Arthropoda"),
    class = c("Copepoda", NA_character_, "Branchiopoda"),
    order = c("Calanoida", NA_character_, NA_character_),
    genus = c(NA_character_, NA_character_, "Daphnia"),
    stringsAsFactors = FALSE
  )
  evidence <- rbind(
    calanoida_like_evidence(),
    data.frame(
      taxon_name = "Oomycota",
      taxon_rank = "phylum",
      source = "literature",
      evidence_type = "ecology",
      evidence_summary = "Oomycota includes host-associated parasitic and pathogenic taxa with marine members.",
      reference = "Oomycota ecology",
      environment = "marine",
      habitat = "host-associated; parasitic",
      region = NA_character_,
      reference_url = NA_character_,
      stringsAsFactors = FALSE
    )
  )

  possible <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "llm_possible_exclusion",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  possible_prompt <- paste(readLines(attr(possible, "prompt_path", exact = TRUE), warn = FALSE), collapse = "\n")
  expect_true(grepl("Oomycota\tphylum", possible_prompt, fixed = TRUE))
  expect_false(grepl("Calanoida\torder", possible_prompt, fixed = TRUE))
  expect_false(grepl("Daphnia\tgenus", possible_prompt, fixed = TRUE))
  expect_equal(unique(possible$taxon_name[possible$llm_prompt_selected]), "Oomycota")

  flagged <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "llm_flagged",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  flagged_prompt <- paste(readLines(attr(flagged, "prompt_path", exact = TRUE), warn = FALSE), collapse = "\n")
  expect_true(grepl("Oomycota\tphylum", flagged_prompt, fixed = TRUE))
  expect_true(grepl("Daphnia\tgenus", flagged_prompt, fixed = TRUE))
  expect_false(grepl("Calanoida\torder", flagged_prompt, fixed = TRUE))

  evidence_only <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  evidence_only_prompt <- paste(readLines(attr(evidence_only, "prompt_path", exact = TRUE), warn = FALSE), collapse = "\n")
  expect_true(grepl("No query taxa were selected for LLM/manual judgement", evidence_only_prompt, fixed = TRUE))
  expect_false(any(evidence_only$llm_prompt_selected))
})

test_that("apply_flag_taxa_llm_result reads JSON, validates query IDs, and preserves unselected rows", {
  testthat::skip_if_not_installed("jsonlite")
  evidence_first <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = valid_taxon_evidence()[1, , drop = FALSE],
    judgement_mode = "llm_missing_evidence",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  selected_id <- unique(evidence_first$query_id[evidence_first$llm_selected])
  json_result <- data.frame(
    query_id = selected_id,
    taxon_name = "Daphnia",
    taxon_rank = "genus",
    expected_environment_status = "unknown",
    expected_habitat_status = "no_known_habitat_evidence",
    expected_region_status = "not_assessed",
    ecological_status = "unknown",
    occurrence_interpretation = "uncertain",
    recommended_action = "flag_for_review",
    rationale = "Manual JSON rationale.",
    references = "Manual JSON reference.",
    evidence_basis = "conservative_reasoning_only",
    llm_tool_used = FALSE,
    llm_tool_name = NA_character_,
    llm_tool_query = NA_character_,
    llm_tool_evidence_summary = NA_character_,
    stringsAsFactors = FALSE
  )
  json_path <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(json_result, dataframe = "rows", auto_unbox = TRUE, na = "null"), json_path)

  final <- apply_flag_taxa_llm_result(evidence_first, llm_json_path = json_path)

  expect_equal(final$rationale[[2]], "Manual JSON rationale.")
  expect_equal(final$recommended_action[[2]], "flag_for_review")
  expect_equal(final$recommended_action[[1]], evidence_first$recommended_action[[1]])
  expect_equal(final$feature_id, evidence_first$feature_id)
})

test_that("apply_flag_taxa_llm_result repairs missing tool references", {
  evidence_first <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = valid_taxon_evidence()[1, , drop = FALSE],
    judgement_mode = "llm_missing_evidence",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  selected_id <- unique(evidence_first$query_id[evidence_first$llm_selected])
  json_result <- data.frame(
    query_id = selected_id,
    taxon_name = "Daphnia",
    taxon_rank = "genus",
    expected_environment_status = "unknown",
    expected_habitat_status = "no_known_habitat_evidence",
    expected_region_status = "not_assessed",
    ecological_status = "unknown",
    occurrence_interpretation = "uncertain",
    recommended_action = "flag_for_review",
    rationale = "Manual JSON tool rationale.",
    references = "Manual JSON tool reference.",
    evidence_basis = "tool_evidence",
    llm_tool_used = TRUE,
    llm_tool_name = "scite",
    llm_tool_query = "Daphnia marine deep sea",
    llm_tool_evidence_summary = "Tool was queried but returned no verifiable references.",
    llm_tool_references = "",
    stringsAsFactors = FALSE
  )

  final <- apply_flag_taxa_llm_result(evidence_first, llm_result = json_result)

  expect_equal(
    final$llm_tool_references[[2]],
    "Tool used but no verifiable references were returned."
  )
  expect_equal(final$llm_tool_name[[2]], "scite")
})

test_that("apply_flag_taxa_llm_result repairs required missing tool attempts", {
  evidence_first <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = valid_taxon_evidence()[1, , drop = FALSE],
    judgement_mode = "llm_missing_evidence",
    allow_llm_tools = TRUE,
    prompt_tools = "Scite",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  selected_id <- unique(evidence_first$query_id[evidence_first$llm_selected])
  json_result <- data.frame(
    query_id = selected_id,
    taxon_name = "Daphnia",
    taxon_rank = "genus",
    expected_environment_status = "unknown",
    expected_habitat_status = "no_known_habitat_evidence",
    expected_region_status = "not_assessed",
    ecological_status = "unknown",
    occurrence_interpretation = "uncertain",
    recommended_action = "exclude",
    rationale = "Manual JSON attempted exclusion without tool trace.",
    references = "Manual JSON reference.",
    evidence_basis = "conservative_reasoning_only",
    llm_tool_used = FALSE,
    llm_tool_name = NA_character_,
    llm_tool_query = NA_character_,
    llm_tool_evidence_summary = NA_character_,
    llm_tool_references = NA_character_,
    stringsAsFactors = FALSE
  )

  final <- apply_flag_taxa_llm_result(evidence_first, llm_result = json_result)

  expect_true(final$llm_tool_used[[2]])
  expect_equal(final$llm_tool_name[[2]], "Scite")
  expect_equal(
    final$llm_tool_references[[2]],
    "Tool evidence was required but no verifiable references were returned."
  )
  expect_equal(final$recommended_action[[2]], "flag_for_review")
})

test_that("apply_flag_taxa_llm_result errors when selected query IDs are missing", {
  evidence_first <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    judgement_mode = "llm_all",
    prompt_only = TRUE,
    prompt_path = tempdir()
  )
  selected_id <- evidence_first$query_id[[1]]
  llm_result <- data.frame(
    query_id = selected_id,
    taxon_name = evidence_first$taxon_name[[1]],
    taxon_rank = evidence_first$taxon_rank[[1]],
    expected_environment_status = "unknown",
    expected_habitat_status = "unknown",
    expected_region_status = "not_assessed",
    ecological_status = "unknown",
    occurrence_interpretation = "uncertain",
    recommended_action = "flag_for_review",
    rationale = "Only one row.",
    references = "Manual reference.",
    evidence_basis = "conservative_reasoning_only",
    llm_tool_used = FALSE,
    llm_tool_name = NA_character_,
    llm_tool_query = NA_character_,
    llm_tool_evidence_summary = NA_character_,
    stringsAsFactors = FALSE
  )

  expect_error(
    apply_flag_taxa_llm_result(evidence_first, llm_result = llm_result),
    "one row per unique query taxon"
  )
})

test_that("final output column order starts with the documented flag_taxa columns", {
  output <- flag_taxa(
    calanoida_like_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = calanoida_like_evidence(),
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expected_prefix <- c(
    "feature_id",
    "taxon_name",
    "taxon_rank",
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
    "evidence_summary",
    "llm_prompt_selected",
    "llm_tool_used",
    "llm_tool_name",
    "llm_tool_query",
    "llm_tool_evidence_summary",
    "llm_tool_references",
    "expected_environment",
    "expected_habitat",
    "expected_region"
  )
  expect_equal(colnames(output)[seq_along(expected_prefix)], expected_prefix)
})

test_that("deterministic retain rows bypass non-all LLM modes", {
  for (mode in c("llm_missing_evidence", "llm_flagged", "llm_possible_exclusion", "llm_review")) {
    testthat::local_mocked_bindings(
      call_flag_taxa_ellmer = function(...) {
        stop("LLM should not be called for deterministic retain rows.", call. = FALSE)
      },
      call_flag_taxa_tool_evidence = function(...) {
        stop("Tools should not be called for deterministic retain rows.", call. = FALSE)
      }
    )

    output <- flag_taxa(
      calanoida_like_taxonomy(),
      expected_environment = "marine",
      expected_habitat = "deep sea",
      expected_region = "Clarion-Clipperton Zone",
      taxon_evidence = calanoida_like_evidence(),
      judgement_mode = mode,
      allow_llm_tools = TRUE,
      prompt_only = FALSE,
      chat = "fake_chat"
    )

    expect_equal(output$recommended_action, "retain")
    expect_equal(output$evidence_sources, "worms; worms_ccz")
  }
})

test_that("evidence_only can exclude species with explicit environment incompatibility", {
  tax <- data.frame(
    feature_id = "asv_fresh",
    phylum = "Chordata",
    genus = "Freshus",
    species = "Freshus strictus",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Freshus strictus",
    taxon_rank = "species",
    source = "literature",
    evidence_type = "environment",
    evidence_summary = "Species is documented from freshwater and terrestrial habitats only.",
    reference = "Freshwater monograph",
    environment = "freshwater; terrestrial",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    taxon_evidence = evidence,
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_equal(output$expected_environment_status, "incompatible")
  expect_equal(output$recommended_action, "exclude")
  expect_equal(output$evidence_basis, "local_evidence")
})

test_that("evidence_only treats missing evidence as review not exclusion", {
  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    judgement_mode = "evidence_only",
    prompt_only = FALSE
  )

  expect_true(all(output$recommended_action == "flag_for_review"))
  expect_false(any(output$recommended_action == "exclude"))
  expect_true(all(output$expected_environment_status == "unknown"))
  expect_true(all(output$expected_habitat_status == "no_known_habitat_evidence"))
  expect_true(all(grepl("Missing evidence", output$rationale, fixed = TRUE)))
  expect_false(any(grepl("positive ecological concern", output$rationale, fixed = TRUE)))
})

test_that("judgement modes select LLM taxa from preliminary evidence", {
  evidence <- valid_taxon_evidence()[1, , drop = FALSE]
  calls <- list()

  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(prompt, ...) {
      calls[[length(calls) + 1L]] <<- prompt
      if (grepl("q0001", prompt, fixed = TRUE) && grepl("q0002", prompt, fixed = TRUE)) {
        return(add_query_ids(
          valid_flag_taxa_taxon_result(
            expected_region_status = rep("not_assessed", 3)
          )[1:2, , drop = FALSE]
        ))
      }
      data.frame(
        query_id = "q0002",
        taxon_name = "Daphnia",
        taxon_rank = "genus",
        expected_environment_status = "unknown",
        expected_habitat_status = "no_known_habitat_evidence",
        expected_region_status = "not_assessed",
        recommended_action = "flag_for_review",
        rationale = "LLM reviewed missing evidence.",
        references = "Tool reference for Daphnia.",
        evidence_basis = "conservative_reasoning_only",
        llm_tool_used = FALSE,
        llm_tool_name = NA_character_,
        llm_tool_query = NA_character_,
        llm_tool_evidence_summary = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  )

  missing_output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence,
    judgement_mode = "llm_missing_evidence",
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(length(calls), 1)
  expect_false(grepl("Salmo salar\t", calls[[1]], fixed = TRUE))
  expect_true(grepl("Daphnia\t", calls[[1]], fixed = TRUE))
  expect_equal(missing_output$rationale[[2]], "LLM reviewed missing evidence.")

  calls <- list()
  all_output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence,
    judgement_mode = "llm_all",
    prompt_only = FALSE,
    chat = "fake_chat"
  )
  expect_equal(length(calls), 1)
  expect_true(grepl("Salmo salar\t", calls[[1]], fixed = TRUE))
  expect_true(grepl("Daphnia\t", calls[[1]], fixed = TRUE))
  expect_equal(nrow(all_output), 2)
})

test_that("tool pass follows judgement_mode selection", {
  evidence <- valid_taxon_evidence()[1, , drop = FALSE]
  tool_taxa <- character()
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(query_table, ...) {
      tool_taxa <<- c(tool_taxa, query_table$query_name)
      flag_taxa_tool_evidence_result(query_table)
    },
    call_flag_taxa_ellmer = function(prompt, ...) {
      data.frame(
        query_id = "q0002",
        taxon_name = "Daphnia",
        taxon_rank = "genus",
        expected_environment_status = "unknown",
        expected_habitat_status = "no_known_habitat_evidence",
        expected_region_status = "not_assessed",
        recommended_action = "flag_for_review",
        rationale = "LLM reviewed missing evidence.",
        references = "Tool reference for Daphnia.",
        evidence_basis = "tool_evidence",
        llm_tool_used = TRUE,
        llm_tool_name = "scite",
        llm_tool_query = "Daphnia marine",
        llm_tool_evidence_summary = "Explicit tool evidence for Daphnia",
        llm_tool_references = "Daphnia tool reference. doi:10.0000/daphnia",
        stringsAsFactors = FALSE
      )
    }
  )

  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    taxon_evidence = evidence,
    judgement_mode = "llm_missing_evidence",
    allow_llm_tools = TRUE,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(tool_taxa, "Daphnia")
  expect_equal(output$llm_tool_name[[2]], "scite")
  expect_true(grepl("Explicit tool evidence for Daphnia", output$llm_tool_evidence_summary[[2]], fixed = TRUE))
  expect_true(grepl("Daphnia tool reference", output$llm_tool_references[[2]], fixed = TRUE))
})

test_that("llm_flagged uses tools for flagged taxa but not deterministic retain taxa", {
  tool_taxa <- character()
  llm_taxa <- character()
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(query_table, ...) {
      tool_taxa <<- c(tool_taxa, query_table$query_name)
      flag_taxa_tool_evidence_result(query_table)
    },
    call_flag_taxa_ellmer = function(prompt, ...) {
      if (grepl("Daphnia\tgenus", prompt, fixed = TRUE)) {
        llm_taxa <<- c(llm_taxa, "Daphnia")
      }
      if (grepl("Calanoida\torder", prompt, fixed = TRUE)) {
        llm_taxa <<- c(llm_taxa, "Calanoida")
      }
      data.frame(
        query_id = "q0002",
        taxon_name = "Daphnia",
        taxon_rank = "genus",
        expected_environment_status = "unknown",
        expected_habitat_status = "no_known_habitat_evidence",
        expected_region_status = "unknown",
        recommended_action = "flag_for_review",
        rationale = "LLM reviewed only the review taxon.",
        references = "Tool reference for Daphnia.",
        evidence_basis = "tool_evidence",
        llm_tool_used = TRUE,
        llm_tool_name = "scite",
        llm_tool_query = "Daphnia marine deep sea CCZ",
        llm_tool_evidence_summary = "Explicit tool evidence for Daphnia",
        llm_tool_references = "Daphnia deep-sea tool reference. doi:10.0000/daphnia-deep",
        stringsAsFactors = FALSE
      )
    }
  )

  output <- flag_taxa(
    calanoida_like_taxonomy(include_review_taxon = TRUE),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = calanoida_like_evidence(),
    judgement_mode = "llm_flagged",
    allow_llm_tools = TRUE,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(tool_taxa, "Daphnia")
  expect_equal(llm_taxa, "Daphnia")
  expect_equal(output$recommended_action[[1]], "retain")
  expect_equal(output$recommended_action[[2]], "flag_for_review")
  expect_equal(output$llm_tool_name[[2]], "scite")
})

test_that("llm_possible_exclusion uses tools only for possible-exclusion taxa", {
  tax <- data.frame(
    feature_id = c("asv_calanoida", "asv_oomycota", "asv_daphnia"),
    phylum = c("Arthropoda", "Oomycota", "Arthropoda"),
    class = c("Copepoda", NA_character_, "Branchiopoda"),
    order = c("Calanoida", NA_character_, NA_character_),
    genus = c(NA_character_, NA_character_, "Daphnia"),
    stringsAsFactors = FALSE
  )
  evidence <- rbind(
    calanoida_like_evidence(),
    data.frame(
      taxon_name = "Oomycota",
      taxon_rank = "phylum",
      source = "literature",
      evidence_type = "ecology",
      evidence_summary = "Oomycota includes host-associated parasitic and pathogenic taxa with marine members.",
      reference = "Oomycota ecology",
      environment = "marine",
      habitat = "host-associated; parasitic",
      region = NA_character_,
      reference_url = NA_character_,
      stringsAsFactors = FALSE
    )
  )
  tool_taxa <- character()
  llm_prompt <- NULL
  testthat::local_mocked_bindings(
    call_flag_taxa_tool_evidence = function(query_table, ...) {
      tool_taxa <<- c(tool_taxa, query_table$query_name)
      flag_taxa_tool_evidence_result(query_table)
    },
    call_flag_taxa_ellmer = function(prompt, ...) {
      llm_prompt <<- prompt
      data.frame(
        query_id = "q0002",
        taxon_name = "Oomycota",
        taxon_rank = "phylum",
        expected_environment_status = "compatible",
        expected_habitat_status = "no_known_habitat_evidence",
        expected_region_status = "no_distribution_evidence",
        ecological_status = "unlikely_resident",
        occurrence_interpretation = "possible_transient_or_allochthonous",
        recommended_action = "flag_possible_exclusion",
        rationale = "Host-associated/parasitic ecology makes free-living resident deep-sea interpretation unlikely.",
        references = "Oomycota ecology; tool reference.",
        evidence_basis = "local_and_tool_evidence",
        llm_tool_used = TRUE,
        llm_tool_name = "scite",
        llm_tool_query = "Oomycota marine deep sea",
        llm_tool_evidence_summary = "Explicit tool evidence for Oomycota",
        llm_tool_references = "Oomycota tool reference. doi:10.0000/oomycota",
        stringsAsFactors = FALSE
      )
    }
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea benthic sediment",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    judgement_mode = "llm_possible_exclusion",
    allow_llm_tools = TRUE,
    prompt_only = FALSE,
    chat = "fake_chat"
  )

  expect_equal(tool_taxa, "Oomycota")
  expect_true(grepl("Oomycota\tphylum", llm_prompt, fixed = TRUE))
  expect_false(grepl("Calanoida\torder", llm_prompt, fixed = TRUE))
  expect_false(grepl("Daphnia\tgenus", llm_prompt, fixed = TRUE))
  expect_equal(output$recommended_action[[1]], "retain")
  expect_equal(output$recommended_action[[2]], "flag_possible_exclusion")
  expect_equal(output$recommended_action[[3]], "flag_for_review")
  expect_false(output$llm_tool_used[[1]])
  expect_true(output$llm_tool_used[[2]])
  expect_false(output$llm_tool_used[[3]])
})

test_that("prompt includes preliminary judgement baseline instructions", {
  prompt <- flag_taxa_prompt_for_test(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    judgement_mode = "llm_all",
    prompt_only = TRUE
  )

  expect_true(grepl("Preliminary evidence-first judgement", prompt, fixed = TRUE))
  expect_true(grepl("deterministic local-evidence baseline", prompt, fixed = TRUE))
  expect_true(grepl("do not override review to exclude", tolower(prompt), fixed = TRUE))
  expect_true(grepl("possible_likely", prompt, fixed = TRUE))
  expect_true(grepl("possible_unlikely", prompt, fixed = TRUE))
  expect_true(grepl("your rationale must explain the specific biological/ecological reason", prompt, fixed = TRUE))
  expect_true(grepl("Do not merely repeat status values", prompt, fixed = TRUE))
  expect_true(grepl("Absence of records alone is not incompatibility", prompt, fixed = TRUE))
  expect_true(grepl("photosynthetic taxa in aphotic deep-sea benthic habitats", prompt, fixed = TRUE))
})

test_that("repair_invalid_llm downgrades unsupported no-evidence exclusions by default", {
  result <- valid_flag_taxa_taxon_result_with_tool_fields()
  result$expected_habitat_status[[2]] <- "incompatible"
  result$recommended_action[[2]] <- "exclude"
  result$evidence_basis[[2]] <- "tool_evidence"
  result$llm_tool_used[[2]] <- TRUE
  result$llm_tool_evidence_summary[[2]] <- "no_evidence_found: no explicit evidence found for deep sea or CCZ."
  result$references[[2]] <- "No explicit tool evidence."

  output <- flag_taxa(
    flag_taxa_duplicate_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    prompt_only = FALSE,
    llm_result = result
  )

  expect_equal(output$recommended_action[[3]], "flag_for_review")
  expect_equal(output$expected_habitat_status[[3]], "no_known_habitat_evidence")
  expect_true(grepl("biohelper repair", output$rationale[[3]], fixed = TRUE))
})

test_that("no evidence tool summaries cannot justify habitat incompatible or exclude", {
  query_table <- biohelper:::.unique_flag_taxa_query_table(
    biohelper:::.prepare_flag_taxa_feature_query_map(
      extract_tax_table(flag_taxa_duplicate_taxonomy()),
      tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
    )
  )
  result <- add_query_ids(valid_flag_taxa_taxon_result_with_tool_fields())
  result$expected_habitat_status[[2]] <- "incompatible"
  result$recommended_action[[2]] <- "flag_for_review"
  result$evidence_basis[[2]] <- "tool_evidence"
  result$llm_tool_used[[2]] <- TRUE
  result$llm_tool_evidence_summary[[2]] <- "no_evidence_found: no explicit evidence found for deep sea or CCZ."
  result$references[[2]] <- "No explicit tool evidence."

  expect_error(
    biohelper:::validate_flag_taxa_taxon_output(
      result = result,
      query_taxonomy = query_table,
      expected_region = "Clarion-Clipperton Zone",
      allow_llm_tools = FALSE
    ),
    "expected_habitat_status = \"incompatible\".*no evidence found"
  )

  result$expected_habitat_status[[2]] <- "no_known_habitat_evidence"
  result$recommended_action[[2]] <- "exclude"
  expect_error(
    biohelper:::validate_flag_taxa_taxon_output(
      result = result,
      query_taxonomy = query_table,
      expected_region = "Clarion-Clipperton Zone",
      allow_llm_tools = FALSE
    ),
    "recommended_action = \"exclude\".*no_known_habitat_evidence"
  )
})

test_that("Arcellinida-like no-evidence tool result remains review, not exclude", {
  tax <- data.frame(
    feature_id = "asv_arcellinida",
    phylum = "Amoebozoa",
    class = "Tubulinea",
    order = "Arcellinida",
    stringsAsFactors = FALSE
  )
  evidence <- data.frame(
    taxon_name = "Arcellinida",
    taxon_rank = "order",
    source = "worms",
    evidence_type = "taxonomic_environment_database",
    evidence_summary = "WoRMS environment flags include marine; brackish; freshwater; terrestrial.",
    reference = "WoRMS",
    accepted_name = "Arcellinida",
    environment = "marine; brackish; freshwater; terrestrial",
    stringsAsFactors = FALSE
  )
  result <- data.frame(
    query_id = "q0001",
    taxon_name = "Arcellinida",
    taxon_rank = "order",
    expected_environment_status = "mixed_within_rank",
    expected_habitat_status = "no_known_habitat_evidence",
    expected_region_status = "no_distribution_evidence",
    recommended_action = "flag_for_review",
    rationale = "Marine is included among mixed environment flags, but no explicit deep-sea/CCZ evidence was found; absence is uncertainty.",
    references = "WoRMS; no_evidence_found tool search",
    evidence_basis = "local_and_tool_evidence",
    llm_tool_used = TRUE,
    llm_tool_name = "scite",
    llm_tool_query = "Arcellinida marine deep sea Clarion-Clipperton Zone",
    llm_tool_evidence_summary = "no_evidence_found: Tool search found no explicit evidence for Arcellinida in the expected habitat/region; this is treated as missing evidence, not incompatibility.",
    llm_tool_references = "Tool used but no verifiable references were returned.",
    stringsAsFactors = FALSE
  )

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    expected_habitat = "deep sea",
    expected_region = "Clarion-Clipperton Zone",
    taxon_evidence = evidence,
    prompt_only = FALSE,
    llm_result = result,
    allow_llm_tools = TRUE
  )

  expect_equal(output$recommended_action, "flag_for_review")
  expect_equal(output$expected_habitat_status, "no_known_habitat_evidence")
  expect_false(output$expected_habitat_status == "incompatible")
})

test_that("call_flag_taxa_ellmer converts fake chat list output to a data.frame", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      expect_true(include_feature_id)
      expect_null(expected_region)
      "fake_schema"
    }
  )
  fake_chat <- list(
    chat_structured = function(prompt, type) {
      expect_equal(prompt, "fake prompt")
      expect_equal(type, "fake_schema")
      list(
        as.list(valid_flag_taxa_structured_result()[1, , drop = FALSE]),
        as.list(valid_flag_taxa_structured_result()[2, , drop = FALSE])
      )
    }
  )

  output <- biohelper:::call_flag_taxa_ellmer(
    prompt = "fake prompt",
    chat = fake_chat,
    include_feature_id = TRUE
  )

  expect_s3_class(output, "data.frame")
  expect_equal(output$feature_id, c("asv1", "asv2"))
})

test_that("ellmer schema helper excludes not_assessed when expected_region is provided", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      expect_true(include_feature_id)
      expect_equal(expected_region, "North Atlantic")
      region_status <- biohelper:::.flag_taxa_allowed_values(
        expected_region = expected_region
      )$expected_region_status
      expect_false("not_assessed" %in% region_status)
      expect_true("unknown" %in% region_status)
      "fake_schema"
    }
  )
  fake_chat <- list(
    chat_structured = function(prompt, type) {
      expect_equal(type, "fake_schema")
      valid_flag_taxa_structured_result()
    }
  )

  output <- biohelper:::call_flag_taxa_ellmer(
    prompt = "fake prompt",
    chat = fake_chat,
    include_feature_id = TRUE,
    expected_region = "North Atlantic",
    max_tries = 1,
    retry_sleep = 0
  )

  expect_s3_class(output, "data.frame")
})

test_that("call_flag_taxa_ellmer retries transient ellmer errors then succeeds", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      expect_true(include_feature_id)
      expect_null(expected_region)
      "fake_schema"
    }
  )
  attempts <- 0
  fake_chat <- list(
    provider = "openai",
    model = "fake-model",
    chat_structured = function(prompt, type) {
      attempts <<- attempts + 1
      expect_equal(prompt, "fake prompt")
      expect_equal(type, "fake_schema")
      if (attempts == 1) {
        stop("HTTP 520: temporary upstream failure", call. = FALSE)
      }
      valid_flag_taxa_structured_result()
    }
  )

  output <- biohelper:::call_flag_taxa_ellmer(
    prompt = "fake prompt",
    chat = fake_chat,
    include_feature_id = TRUE,
    max_tries = 3,
    retry_sleep = 0
  )

  expect_equal(attempts, 2)
  expect_s3_class(output, "data.frame")
  expect_equal(output$feature_id, c("asv1", "asv2"))
})

test_that("call_flag_taxa_ellmer reports repeated transient failures clearly", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      expect_true(include_feature_id)
      expect_null(expected_region)
      "fake_schema"
    }
  )
  attempts <- 0
  fake_chat <- list(
    provider = "openai",
    model = "fake-model",
    chat_structured = function(prompt, type) {
      attempts <<- attempts + 1
      stop("HTTP 520: temporary upstream failure", call. = FALSE)
    }
  )

  error <- tryCatch(
    biohelper:::call_flag_taxa_ellmer(
      prompt = "fake prompt",
      chat = fake_chat,
      include_feature_id = TRUE,
      max_tries = 3,
      retry_sleep = 0
    ),
    error = identity
  )

  expect_s3_class(error, "error")
  expect_equal(attempts, 3)
  message <- conditionMessage(error)
  expect_match(message, "failed after 3 attempt")
  expect_match(message, "provider=openai")
  expect_match(message, "model=fake-model")
  expect_match(message, "HTTP 520")
  expect_match(message, "Retry later")
  expect_match(message, "OpenAI status")
  expect_match(message, "chat\\$chat\\(\"hello\"\\)")
})

test_that("HTTP 413 ellmer errors describe payload size rather than quota", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      "fake_schema"
    }
  )
  attempts <- 0
  fake_chat <- list(
    provider = "github",
    model = "gpt-4o-mini",
    chat_structured = function(prompt, type) {
      attempts <<- attempts + 1
      stop("HTTP 413: Payload Too Large", call. = FALSE)
    }
  )

  error <- tryCatch(
    biohelper:::call_flag_taxa_ellmer(
      prompt = "fake prompt",
      chat = fake_chat,
      include_feature_id = FALSE,
      max_tries = 3,
      retry_sleep = 0
    ),
    error = identity
  )

  expect_s3_class(error, "error")
  expect_equal(attempts, 1)
  message <- conditionMessage(error)
  expect_match(message, "HTTP 413")
  expect_match(message, "payload/model context limit")
  expect_match(message, "not provider quota")
  expect_match(message, "prompt_only = TRUE")
  expect_false(grepl("rate-limit", message, fixed = TRUE))
})

test_that("HTTP 429 ellmer errors describe provider quota and not Scite quota", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      "fake_schema"
    }
  )
  fake_chat <- list(
    provider = "github",
    model = "gpt-4o-mini",
    chat_structured = function(prompt, type) {
      stop("HTTP 429: rate limit exceeded", call. = FALSE)
    }
  )

  error <- tryCatch(
    biohelper:::call_flag_taxa_ellmer(
      prompt = "fake prompt",
      chat = fake_chat,
      include_feature_id = FALSE,
      max_tries = 1,
      retry_sleep = 0
    ),
    error = identity
  )

  expect_s3_class(error, "error")
  message <- conditionMessage(error)
  expect_match(message, "HTTP 429")
  expect_match(message, "LLM provider request quota/rate-limit")
  expect_match(message, "does not necessarily mean Scite")
  expect_match(message, "does not necessarily indicate prompt-size problems")
  expect_false(grepl("Scite quota", message, fixed = TRUE))
  expect_false(grepl("payload/model context limit", message, fixed = TRUE))
})

test_that("premature EOF structured-output errors get prompt-size guidance", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() NULL,
    .flag_taxa_ellmer_output_type = function(include_feature_id, expected_region) {
      "fake_schema"
    }
  )
  attempts <- 0
  fake_chat <- list(
    provider = "github",
    model = "gpt-4o-mini",
    chat_structured = function(prompt, type) {
      attempts <<- attempts + 1
      stop("parse error: premature EOF [ { \"query_id\": \"q0001\",", call. = FALSE)
    }
  )

  error <- tryCatch(
    biohelper:::call_flag_taxa_ellmer(
      prompt = "fake prompt",
      chat = fake_chat,
      include_feature_id = FALSE,
      max_tries = 3,
      retry_sleep = 0
    ),
    error = identity
  )

  expect_s3_class(error, "error")
  expect_equal(attempts, 1)
  message <- conditionMessage(error)
  expect_match(message, "truncated or incomplete")
  expect_match(message, "max_evidence_rows_per_query")
  expect_match(message, "larger-context / larger-output model")
})

test_that("large one-chunk final prompt warns before LLM result processing", {
  expect_warning(
    biohelper:::.warn_flag_taxa_large_prompt_chunks(
      prompt_chunks = list(chunk_001 = paste(rep("x", 50), collapse = "")),
      query_chunks = list(chunk_001 = data.frame(query_id = "q0001")),
      warn_prompt_chars = 10,
      llm_chunk_size = Inf
    ),
    "final judgement prompt is [0-9]+ chars in one chunk"
  )
})

test_that("chat backend errors clearly when ellmer is unavailable", {
  testthat::local_mocked_bindings(
    .require_ellmer_for_flag_taxa = function() {
      stop(
        "The ellmer package is required for LLM-backed flag_taxa(). Install it with install.packages('ellmer') or use prompt_only = TRUE.",
        call. = FALSE
      )
    }
  )

  expect_error(
    biohelper:::call_flag_taxa_ellmer(
      prompt = "fake prompt",
      chat = list(chat_structured = function(...) data.frame()),
      include_feature_id = TRUE
    ),
    "ellmer package is required.*install.packages\\('ellmer'\\).*prompt_only = TRUE"
  )
})

test_that("invalid LLM output fails validation", {
  result <- valid_flag_taxa_structured_result()
  result$recommended_action[1] <- "maybe"
  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(...) result
  )

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      chat = "fake_chat"
    ),
    "recommended_action.*invalid values.*maybe"
  )
})

test_that("LLM output with missing required columns fails validation", {
  result <- valid_flag_taxa_structured_result()
  result$references <- NULL
  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(...) result
  )

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      chat = "fake_chat"
    ),
    "missing required flag_taxa columns: references"
  )
})

test_that("LLM output with invalid status values fails validation", {
  result <- valid_flag_taxa_structured_result()
  result$expected_environment_status[1] <- "probably"
  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(...) result
  )

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      chat = "fake_chat"
    ),
    "expected_environment_status.*invalid values.*probably"
  )
})

test_that("LLM output rejects not_assessed region status when expected_region is provided", {
  result <- valid_flag_taxa_structured_result()
  result$expected_region_status[1] <- "not_assessed"
  testthat::local_mocked_bindings(
    call_flag_taxa_ellmer = function(...) result
  )

  expect_error(
    flag_taxa(
      flag_taxa_test_taxonomy(),
      expected_environment = "marine",
      expected_region = "North Atlantic",
      prompt_only = FALSE,
      chat = "fake_chat"
    ),
    "expected_region_status.*must not be `not_assessed`"
  )
})

test_that("mock_llm_result validation does not fetch WoRMS evidence", {
  testthat::local_mocked_bindings(
    .flag_taxa_extract_taxa_for_evidence = function(...) {
      stop("extract_taxa_for_evidence should not be called.", call. = FALSE)
    },
    .flag_taxa_fetch_worms_evidence = function(...) {
      stop("fetch_worms_evidence should not be called.", call. = FALSE)
    }
  )

  result <- valid_flag_taxa_result()
  output <- flag_taxa(
    flag_taxa_test_taxonomy(),
    expected_environment = "marine",
    expected_habitat = "estuary",
    expected_region = "North Atlantic",
    evidence_sources = "worms",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_equal(output$feature_id, result$feature_id)
  expect_equal(output$taxon_name, c("Salmo salar", "Daphnia"))
  expect_equal(output$taxon_rank, c("species", "genus"))
  expect_equal(output$expected_environment, rep("marine", nrow(output)))
  expect_equal(output$expected_habitat, rep("estuary", nrow(output)))
  expect_equal(output$expected_region, rep("North Atlantic", nrow(output)))
  expect_equal(output$evidence_sources, rep("none", nrow(output)))
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

test_that("mock_llm_result with not_assessed passes when expected_region is NULL", {
  tax <- flag_taxa_test_taxonomy()
  result <- valid_flag_taxa_structured_result()
  result$expected_region_status <- "not_assessed"

  output <- flag_taxa(
    tax,
    expected_environment = "marine",
    prompt_only = FALSE,
    mock_llm_result = result
  )

  expect_equal(output$feature_id, result$feature_id)
  expect_equal(output$expected_region_status, rep("not_assessed", nrow(output)))
  expect_true(all(is.na(output$expected_region)))
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

  prompt <- flag_taxa_prompt_for_test(ps, expected_environment = "marine")

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
