test_that("ambiguous NCBI names can be resolved using phylum", {
  input <- data.frame(
    feature_id = c("asv_1", "asv_2"),
    domain = "Eukaryota",
    phylum = c("TargetPhylum", "OtherPhylum"),
    class = "OldClass",
    order = "OldOrder",
    family = "OldFamily",
    genus = "SharedName"
  )
  reference <- data.frame(
    domain = "Eukaryota",
    phylum = c("OtherPhylum", "ThirdPhylum", "TargetPhylum"),
    class = paste0("NewClass", 1:3),
    order = paste0("NewOrder", 1:3),
    family = paste0("NewFamily", 1:3),
    genus = "SharedName",
    row.names = c("11", "22", "33")
  )
  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) {
      rep("11,22,33", length(taxa))
    },
    getRawTaxonomy = function(ids, sqlFile) {
      stats::setNames(lapply(ids, function(id) rev(unlist(reference[id, ]))), ids)
    },
    getTaxonomy = function(ids, sqlFile, desiredTaxa) {
      reference[as.character(ids), desiredTaxa, drop = FALSE]
    },
    .package = "taxonomizr"
  )

  output <- taxo_normalisation(
    input, sqlFile = "mock.sql", addExtra = FALSE, ranks = names(reference)
  )
  expect_equal(output$feature_id, input$feature_id)
  expect_equal(output$phylum, input$phylum)
  expect_equal(output$class, c("NewClass3", "NewClass1"))
  expect_equal(output$genus, input$genus)

  # Phylum must still be used to resolve a name when it is not requested in output.
  output <- taxo_normalisation(
    input, sqlFile = "mock.sql", addExtra = FALSE, ranks = c("class", "genus")
  )
  expect_named(output, c("feature_id", "class", "genus"))
  expect_equal(output$class, c("NewClass3", "NewClass1"))

  # A single feature and missing intermediate ranks must also reach phylum.
  sparse <- input[1, c("feature_id", "domain", "phylum", "genus")]
  output <- taxo_normalisation(
    sparse, sqlFile = "mock.sql", addExtra = FALSE, ranks = names(reference)
  )
  expect_equal(output$phylum, "TargetPhylum")
  expect_equal(output$class, "NewClass3")

  counts <- matrix(1, dimnames = list("asv_1", "sample1"))
  taxonomy <- as.matrix(sparse[, -1, drop = FALSE])
  rownames(taxonomy) <- sparse$feature_id
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy)
  )
  original <- ps
  output <- taxo_normalisation(
    ps, sqlFile = "mock.sql", addExtra = FALSE, ranks = names(reference)
  )
  expect_s4_class(output, "phyloseq")
  expect_equal(unname(phyloseq::tax_table(output)[, "class"]@.Data[1, 1]), "NewClass3")
  expect_identical(ps, original)
})

test_that("all broader names can resolve ambiguity and unsupported names stay unassigned", {
  input <- data.frame(
    feature_id = c("ambiguous", "unmatched", "unique", "absent"),
    domain = "Eukaryota",
    kingdom = c("SharedKingdom", "OtherKingdom", "OtherKingdom", "OtherKingdom"),
    phylum = c("SharedPhylum", "UnmatchedPhylum", "SharedPhylum", "SharedPhylum"),
    genus = c("SharedGenus", "SharedGenus", "UniqueGenus", "AbsentGenus")
  )
  reference <- data.frame(
    domain = "Eukaryota",
    kingdom = "SharedKingdom",
    phylum = "SharedPhylum",
    genus = c("SharedGenus", "SharedGenus"),
    row.names = c("11", "22")
  )
  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) {
      unname(c(SharedGenus = "22,11", UniqueGenus = "11")[taxa])
    },
    getRawTaxonomy = function(ids, sqlFile) {
      stats::setNames(lapply(ids, function(id) rev(unlist(reference[id, ]))), ids)
    },
    getTaxonomy = function(ids, sqlFile, desiredTaxa) {
      expect_false(any(grepl(",", ids, fixed = TRUE), na.rm = TRUE))
      reference[as.character(ids), desiredTaxa, drop = FALSE]
    },
    .package = "taxonomizr"
  )
  expect_warning(
    output <- taxo_normalisation(
      input, sqlFile = "mock.sql", addExtra = FALSE, ranks = names(reference)
    ),
    "Multiple NCBI IDs.*ambiguous.*first matching ID"
  )
  expect_equal(output$kingdom[1], "SharedKingdom")
  # The domain is still available as evidence when closer names are unknown.
  expect_equal(output$genus[2], "SharedGenus")
  expect_true(is.na(output$genus[4]))
  expect_equal(output$domain[4], "Unknown")
  expect_equal(output$genus[3], "SharedGenus")

  # Without a supplied higher rank there is no evidence to choose an ID.
  output <- taxo_normalisation(
    input[1, c("feature_id", "genus")], sqlFile = "mock.sql",
    addExtra = FALSE, ranks = "genus"
  )
  expect_true(is.na(output$genus))
})

test_that("a unique match at a closer rank resolves a name before phylum", {
  input <- data.frame(
    feature_id = "asv_1", phylum = "OldPhylum", family = "TargetFamily", genus = "SharedGenus"
  )
  reference <- data.frame(
    phylum = "NewPhylum", family = c("OtherFamily", "TargetFamily"),
    genus = "SharedGenus", row.names = c("11", "22")
  )
  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) "11,22",
    getRawTaxonomy = function(ids, sqlFile) {
      stats::setNames(lapply(ids, function(id) rev(unlist(reference[id, ]))), ids)
    },
    getTaxonomy = function(ids, sqlFile, desiredTaxa) {
      reference[as.character(ids), desiredTaxa, drop = FALSE]
    },
    .package = "taxonomizr"
  )
  output <- taxo_normalisation(
    input, sqlFile = "mock.sql", addExtra = FALSE, ranks = names(reference)
  )
  expect_equal(output$family, "TargetFamily")
  expect_equal(output$phylum, "NewPhylum")
})

test_that("multiple matching candidates are retained for the next broader clue", {
  input <- data.frame(
    feature_id = "asv_1", phylum = "TargetPhylum", family = "SharedGroup", genus = "SharedName"
  )
  lineages <- list(
    `11` = c(genus = "SharedName", superfamily = "SharedGroup", phylum = "OtherPhylum"),
    `22` = c(genus = "SharedName", superfamily = "SharedGroup", phylum = "TargetPhylum"),
    `33` = c(genus = "SharedName", superfamily = "OtherGroup", phylum = "TargetPhylum")
  )
  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) "11,22,33",
    getRawTaxonomy = function(ids, sqlFile) lineages[ids],
    getTaxonomy = function(ids, sqlFile, desiredTaxa) {
      expect_equal(ids, "22")
      matrix(c("TargetPhylum", "SharedName"), nrow = 1,
             dimnames = list("22", desiredTaxa))
    },
    .package = "taxonomizr"
  )
  expect_warning(
    output <- taxo_normalisation(
      input, sqlFile = "mock.sql", addExtra = FALSE, ranks = c("phylum", "genus")
    ),
    regexp = NA
  )
  expect_equal(output$phylum, "TargetPhylum")
  expect_equal(output$genus, "SharedName")
})
