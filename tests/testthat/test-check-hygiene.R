test_that("DNAStringSet_to_df works for DNAStringSet input", {
  dss <- Biostrings::DNAStringSet(c(seq1 = "ACGT", seq2 = "AC"))

  output <- DNAStringSet_to_df(dss)

  expect_equal(output$width, c(4L, 2L))
  expect_equal(output$seq, c("ACGT", "AC"))
  expect_equal(output$names, c("seq1", "seq2"))
})

test_that("taxo_normalisation does not warn about deprecated across arguments", {
  tax_table <- data.frame(
    feature_id = "asv_1",
    superkingdom = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    genus = "Homo",
    species = "Homo_sapiens",
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) {
      expect_equal(taxa, "Homo sapiens")
      "9606"
    },
    getTaxonomy = function(ids, sqlFile, desiredTaxa) {
      expect_equal(ids, "9606")
      stats::setNames(
        data.frame(
          Eukaryota = "Eukaryota",
          Metazoa = "Metazoa",
          Chordata = "Chordata",
          Homo = "Homo",
          `Homo sapiens` = "Homo sapiens",
          check.names = FALSE
        ),
        desiredTaxa
      )
    },
    .package = "taxonomizr"
  )

  expect_warning(
    output <- taxo_normalisation(
      tax_table,
      sqlFile = "mock.sql",
      ranks = c("Superkingdom", "Kingdom", "Phylum", "Genus", "Species"),
      addExtra = FALSE
    ),
    regexp = NA
  )

  expect_equal(output$species, "Homo sapiens")
})

test_that("tax_fix handles unknowns = NULL", {
  tax_table <- data.frame(
    Phylum = c("Chordata", "unknown"),
    Family = c("f__", "Felidae"),
    Genus = c("Panthera", "unknown"),
    row.names = c("taxon1", "taxon2"),
    stringsAsFactors = FALSE
  )

  output <- tax_fix(tax_table, unknowns = NULL, verbose = FALSE)

  expect_s3_class(output, "data.frame")
  expect_equal(output$Genus[2], "unknown")
})

test_that("pw_group_dissimilarity works without external helper functions", {
  otu <- matrix(
    c(
      1, 0, 1, 0,
      0, 1, 0, 1,
      1, 1, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(otu) <- paste0("taxon", 1:3)
  colnames(otu) <- paste0("sample", 1:4)
  metadata <- data.frame(
    Group = c("A", "A", "B", "B"),
    row.names = colnames(otu)
  )
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(metadata)
  )

  output <- pw_group_dissimilarity(ps, group = "Group", method = "jaccard", justDF = TRUE)

  expect_s3_class(output, "data.frame")
  expect_named(output, c("Group", "Dist", "Comparison"))
  expect_true("between-groups" %in% output$Comparison)
})
