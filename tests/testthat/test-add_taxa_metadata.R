add_taxa_metadata_test_ps <- function() {
  otu <- phyloseq::otu_table(
    matrix(
      c(10, 0, 2, 3, 1, 8),
      nrow = 3,
      dimnames = list(c("asv1", "asv2", "asv3"), c("sample1", "sample2"))
    ),
    taxa_are_rows = TRUE
  )
  tax <- phyloseq::tax_table(
    matrix(
      c(
        "Metazoa", "Chordata", "Salmo",
        "Metazoa", "Arthropoda", "Daphnia",
        "Fungi", "Ascomycota", "Candida"
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("asv1", "asv2", "asv3"), c("kingdom", "phylum", "genus"))
    )
  )
  phyloseq::phyloseq(otu, tax)
}

add_taxa_metadata_test_list <- function() {
  list(
    summary = data.frame(
      feature_id = c("asv2", "asv1"),
      taxon_name = c("Daphnia", "Salmo"),
      recommended_action = c("flag_for_review", "retain"),
      rationale = c("needs review", "marine evidence"),
      expected_environment_status = c("incompatible", "compatible"),
      stringsAsFactors = FALSE
    ),
    obis = data.frame(
      feature_id = c("asv1", "asv3"),
      obis_region_status = c("obis_records_found", "obis_not_queried"),
      obis_total_records_including_descendants = c(12L, 0L),
      stringsAsFactors = FALSE
    )
  )
}

test_that("add_taxa_metadata default adds all summary columns except feature_id", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    verbose = FALSE
  )
  tax <- as(phyloseq::tax_table(ps), "matrix")

  expect_true(all(c(
    "tm_taxon_name",
    "tm_recommended_action",
    "tm_rationale",
    "tm_expected_environment_status"
  ) %in% colnames(tax)))
  expect_false("tm_feature_id" %in% colnames(tax))
})

test_that("add_taxa_metadata prefixes added columns with tm_", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )

  expect_true("tm_recommended_action" %in% phyloseq::rank_names(ps))
})

test_that("add_taxa_metadata accepts selected cols", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = c("recommended_action", "rationale"),
    verbose = FALSE
  )
  tax <- as(phyloseq::tax_table(ps), "matrix")

  expect_true(all(c("tm_recommended_action", "tm_rationale") %in% colnames(tax)))
  expect_false("tm_taxon_name" %in% colnames(tax))
})

test_that("add_taxa_metadata can select another source table", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    source = "obis",
    cols = "obis_region_status",
    verbose = FALSE
  )
  tax <- as(phyloseq::tax_table(ps), "matrix")

  expect_equal(tax["asv1", "tm_obis_region_status"], "obis_records_found")
  expect_equal(tax["asv3", "tm_obis_region_status"], "obis_not_queried")
})

test_that("add_taxa_metadata accepts a full assess_taxa_evidence-style list", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )

  expect_equal(
    as(phyloseq::tax_table(ps), "matrix")["asv1", "tm_recommended_action"],
    "retain"
  )
})

test_that("add_taxa_metadata accepts a single metadata data frame", {
  metadata <- add_taxa_metadata_test_list()$summary
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    metadata,
    cols = "recommended_action",
    verbose = FALSE
  )

  expect_equal(
    as(phyloseq::tax_table(ps), "matrix")["asv2", "tm_recommended_action"],
    "flag_for_review"
  )
})

test_that("add_taxa_metadata matches by feature_id against taxa_names", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )
  tax <- as(phyloseq::tax_table(ps), "matrix")

  expect_equal(tax["asv1", "tm_recommended_action"], "retain")
  expect_equal(tax["asv2", "tm_recommended_action"], "flag_for_review")
})

test_that("add_taxa_metadata preserves taxa and inserts NA for missing metadata rows", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )
  tax <- as(phyloseq::tax_table(ps), "matrix")

  expect_equal(phyloseq::taxa_names(ps), c("asv1", "asv2", "asv3"))
  expect_true(is.na(tax["asv3", "tm_recommended_action"]))
})

test_that("add_taxa_metadata errors on duplicate feature_id values", {
  metadata <- data.frame(
    feature_id = c("asv1", "asv1"),
    recommended_action = c("retain", "exclude"),
    stringsAsFactors = FALSE
  )

  expect_error(
    add_taxa_metadata(add_taxa_metadata_test_ps(), metadata, verbose = FALSE),
    "Duplicate `feature_id`"
  )
})

test_that("add_taxa_metadata errors on missing requested columns", {
  expect_error(
    add_taxa_metadata(
      add_taxa_metadata_test_ps(),
      add_taxa_metadata_test_list(),
      cols = c("recommended_action", "missing_column"),
      verbose = FALSE
    ),
    "not found"
  )
})

test_that("add_taxa_metadata errors on existing metadata columns by default", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )

  expect_error(
    add_taxa_metadata(ps, add_taxa_metadata_test_list(), cols = "recommended_action", verbose = FALSE),
    "already exist"
  )
})

test_that("add_taxa_metadata overwrites existing metadata columns when requested", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )
  replacement <- data.frame(
    feature_id = c("asv1", "asv2"),
    recommended_action = c("exclude", "retain"),
    stringsAsFactors = FALSE
  )
  ps <- add_taxa_metadata(
    ps,
    replacement,
    cols = "recommended_action",
    overwrite = TRUE,
    verbose = FALSE
  )

  expect_equal(
    as(phyloseq::tax_table(ps), "matrix")["asv1", "tm_recommended_action"],
    "exclude"
  )
})

test_that("add_taxa_metadata preserves existing taxonomy columns exactly", {
  ps <- add_taxa_metadata_test_ps()
  before <- as(phyloseq::tax_table(ps), "matrix")
  out <- add_taxa_metadata(ps, add_taxa_metadata_test_list(), verbose = FALSE)
  after <- as(phyloseq::tax_table(out), "matrix")

  expect_identical(after[, colnames(before), drop = FALSE], before)
})

test_that("add_taxa_metadata appends metadata columns after taxonomy columns", {
  ps <- add_taxa_metadata_test_ps()
  before_names <- colnames(as(phyloseq::tax_table(ps), "matrix"))
  out <- add_taxa_metadata(
    ps,
    add_taxa_metadata_test_list(),
    cols = c("recommended_action", "rationale"),
    verbose = FALSE
  )
  after_names <- colnames(as(phyloseq::tax_table(out), "matrix"))

  expect_identical(after_names[seq_along(before_names)], before_names)
  expect_identical(tail(after_names, 2), c("tm_recommended_action", "tm_rationale"))
})

test_that("add_taxa_metadata safely converts list columns to character", {
  metadata <- data.frame(
    feature_id = c("asv1", "asv2", "asv3"),
    stringsAsFactors = FALSE
  )
  metadata$evidence <- I(list(c("ref1", "ref2"), "ref3", character()))
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    metadata,
    cols = "evidence",
    verbose = FALSE
  )
  tax <- as(phyloseq::tax_table(ps), "matrix")

  expect_equal(tax["asv1", "tm_evidence"], "ref1; ref2")
  expect_equal(tax["asv2", "tm_evidence"], "ref3")
  expect_true(is.na(tax["asv3", "tm_evidence"]))
})

test_that("drop_taxa_metadata removes columns starting with tm_", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )
  ps <- drop_taxa_metadata(ps)

  expect_false("tm_recommended_action" %in% phyloseq::rank_names(ps))
})

test_that("drop_taxa_metadata retains normal taxonomy columns", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    verbose = FALSE
  )
  before <- as(phyloseq::tax_table(add_taxa_metadata_test_ps()), "matrix")
  ps <- drop_taxa_metadata(ps)

  expect_identical(as(phyloseq::tax_table(ps), "matrix"), before)
})

test_that("drop_taxa_metadata is harmless when no metadata columns are present", {
  ps <- add_taxa_metadata_test_ps()
  out <- drop_taxa_metadata(ps)

  expect_identical(as(phyloseq::tax_table(out), "matrix"), as(phyloseq::tax_table(ps), "matrix"))
})

test_that("drop_taxa_metadata supports a custom prefix", {
  ps <- add_taxa_metadata(
    add_taxa_metadata_test_ps(),
    add_taxa_metadata_test_list(),
    cols = "recommended_action",
    prefix = "taxmeta_",
    verbose = FALSE
  )
  ps <- drop_taxa_metadata(ps, prefix = "taxmeta_")

  expect_false("taxmeta_recommended_action" %in% phyloseq::rank_names(ps))
  expect_true(all(c("kingdom", "phylum", "genus") %in% phyloseq::rank_names(ps)))
})
