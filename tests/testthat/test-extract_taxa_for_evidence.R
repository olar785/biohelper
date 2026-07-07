evidence_taxonomy_table <- function() {
  tax <- data.frame(
    Domain = rep("Eukaryota", 8),
    Superkingdom = rep("Eukaryota", 8),
    Kingdom = c(
      "Animalia", "Animalia", "Bacteria", "Animalia",
      "Animalia", "Bacteria", "Animalia", "Animalia"
    ),
    PHYLUM = c(
      "Chordata", "Chordata", "Proteobacteria", "Chordata",
      "Cnidaria", "", "Chordata", "Chordata"
    ),
    Class = c(
      "Actinopteri", "Actinopteri", "Gammaproteobacteria", "Mammalia",
      "", "", "Actinopteri", "Actinopteri"
    ),
    order = c(
      "Salmoniformes", "Gadiformes", "uncultured", "Cetacea",
      "", "", "Salmoniformes", "Gadiformes"
    ),
    Family = c(
      "Salmonidae", "Gadidae", "Incertae sedis", "Delphinidae",
      "", "", "Salmonidae", "Gadidae"
    ),
    Genus = c("Salmo", "Gadus", "g__", "", "", "", "Salmo", "Gadus"),
    Species = c("Salmo salar", "", NA, "unidentified", "", "", NA, NA),
    stringsAsFactors = FALSE
  )
  rownames(tax) <- paste0("asv", seq_len(nrow(tax)))
  tax
}

load_ps_test_data_euk <- function() {
  testthat::skip_if_not_installed("phyloseq")
  data("ps_test_data_euk", package = "biohelper", envir = environment())
  ps_test_data_euk
}

test_that("mode highest works with a data.frame taxonomy table", {
  taxa <- extract_taxa_for_evidence(evidence_taxonomy_table())

  expect_s3_class(taxa, "tbl_df")
  expect_equal(colnames(taxa), c("taxon_name", "taxon_rank"))
  expect_equal(nrow(taxa), 6)
  expect_equal(anyDuplicated(paste(taxa$taxon_name, taxa$taxon_rank, sep = "\r")), 0)
})

test_that("mode highest skips requested ranks that are absent from the table", {
  tax_table <- data.frame(
    Kingdom = "Animalia",
    Phylum = "Chordata",
    Class = "Actinopteri",
    Family = "Salmonidae",
    Genus = "Salmo",
    Species = "Salmo salar",
    stringsAsFactors = FALSE
  )

  taxa <- extract_taxa_for_evidence(tax_table)

  expect_equal(nrow(taxa), 1)
  expect_equal(taxa$taxon_name, "Salmo salar")
  expect_equal(taxa$taxon_rank, "Species")
})

test_that("mode highest works with a matrix taxonomy table", {
  taxa <- extract_taxa_for_evidence(as.matrix(evidence_taxonomy_table()))

  expect_s3_class(taxa, "tbl_df")
  expect_equal(nrow(taxa), 6)
  expect_true(any(taxa$taxon_name == "Delphinidae" & taxa$taxon_rank == "Family"))
})

test_that("mode highest works with a phyloseq object", {
  testthat::skip_if_not_installed("phyloseq")

  otu <- phyloseq::otu_table(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("asv1", "asv2"), c("sample1", "sample2"))
    ),
    taxa_are_rows = TRUE
  )
  tax <- phyloseq::tax_table(
    matrix(
      c("Animalia", "Chordata", "Salmonidae", "Salmo", "Animalia", "Chordata", "Gadidae", "Gadus"),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("asv1", "asv2"), c("Kingdom", "Phylum", "Family", "Genus"))
    )
  )
  ps <- phyloseq::phyloseq(otu, tax)

  taxa <- extract_taxa_for_evidence(
    ps,
    ranks = c("Family", "Genus"),
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_equal(taxa$feature_id, c("asv1", "asv2"))
  expect_equal(taxa$taxon_name, c("Salmo", "Gadus"))
  expect_equal(taxa$taxon_rank, c("Genus", "Genus"))
})

test_that("mode highest chooses Species when Species is available", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    include_feature_id = TRUE,
    unique = FALSE
  )
  asv1 <- taxa[taxa$feature_id == "asv1", ]

  expect_equal(asv1$taxon_name, "Salmo salar")
  expect_equal(asv1$taxon_rank, "Species")
})

test_that("mode highest falls back to Genus when Species is missing", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    include_feature_id = TRUE,
    unique = FALSE
  )
  asv2 <- taxa[taxa$feature_id == "asv2", ]

  expect_equal(asv2$taxon_name, "Gadus")
  expect_equal(asv2$taxon_rank, "Genus")
})

test_that("mode highest falls back to Family, Class, and Phylum", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    include_feature_id = TRUE,
    unique = FALSE
  )

  asv3 <- taxa[taxa$feature_id == "asv3", ]
  asv4 <- taxa[taxa$feature_id == "asv4", ]
  asv5 <- taxa[taxa$feature_id == "asv5", ]

  expect_equal(asv3$taxon_name, "Gammaproteobacteria")
  expect_equal(asv3$taxon_rank, "Class")
  expect_equal(asv4$taxon_name, "Delphinidae")
  expect_equal(asv4$taxon_rank, "Family")
  expect_equal(asv5$taxon_name, "Cnidaria")
  expect_equal(asv5$taxon_rank, "Phylum")
})

test_that("Domain, Superkingdom, and Kingdom are ignored by default", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_false(any(taxa$taxon_rank %in% c("Domain", "Superkingdom", "Kingdom")))
  expect_false("asv6" %in% taxa$feature_id)
})

test_that("Kingdom can be queried when explicitly requested", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Kingdom",
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_true("asv6" %in% taxa$feature_id)
  expect_true(any(taxa$taxon_name == "Bacteria" & taxa$taxon_rank == "Kingdom"))
})

test_that("Phylum is retained as the broadest default fallback", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_true(any(taxa$feature_id == "asv5" & taxa$taxon_name == "Cnidaria"))
  expect_true(any(taxa$feature_id == "asv5" & taxa$taxon_rank == "Phylum"))
})

test_that("mode all returns multiple rank rows per feature", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    mode = "all",
    include_feature_id = TRUE,
    unique = FALSE
  )
  asv1 <- taxa[taxa$feature_id == "asv1", ]

  expect_true(nrow(asv1) > 1)
  expect_equal(
    asv1$taxon_rank,
    c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
})

test_that("placeholder and empty taxa are dropped", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    mode = "all",
    unique = FALSE
  )

  expect_false(any(is.na(taxa$taxon_name)))
  expect_false(any(taxa$taxon_name %in% c("", "uncultured", "g__", "Incertae sedis", "unidentified")))
})

test_that("include_feature_id controls feature ID output", {
  taxa_with_ids <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Genus",
    include_feature_id = TRUE,
    unique = FALSE
  )
  taxa_without_ids <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Genus",
    include_feature_id = FALSE
  )

  expect_true("feature_id" %in% colnames(taxa_with_ids))
  expect_false("feature_id" %in% colnames(taxa_without_ids))
})

test_that("include_feature_id uses feature ID columns before row names", {
  tax <- data.frame(
    feature_id = c("feature_a", "feature_b"),
    Genus = c("Salmo", "Gadus"),
    stringsAsFactors = FALSE
  )

  taxa <- extract_taxa_for_evidence(tax, ranks = "Genus")

  expect_false("feature_id" %in% colnames(taxa))

  taxa_with_ids <- extract_taxa_for_evidence(
    tax,
    ranks = "Genus",
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_equal(taxa_with_ids$feature_id, c("feature_a", "feature_b"))
})

test_that("unique behaviour is respected in mode highest", {
  taxa_not_unique <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Genus",
    include_feature_id = TRUE,
    unique = FALSE
  )
  taxa_unique <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Genus",
    include_feature_id = TRUE,
    unique = TRUE
  )

  expect_equal(sum(taxa_not_unique$taxon_name == "Gadus"), 2)
  expect_equal(sum(taxa_unique$taxon_name == "Gadus"), 1)
  expect_true("feature_id" %in% colnames(taxa_unique))
})

test_that("unique behaviour is respected in mode all", {
  taxa_not_unique <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Family",
    mode = "all",
    unique = FALSE
  )
  taxa_unique <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Family",
    mode = "all",
    unique = TRUE
  )

  expect_equal(sum(taxa_not_unique$taxon_name == "Gadidae"), 2)
  expect_equal(sum(taxa_unique$taxon_name == "Gadidae"), 1)
})

test_that("unique output ignores feature_id in the deduplication key", {
  taxa <- extract_taxa_for_evidence(
    evidence_taxonomy_table(),
    ranks = "Genus",
    include_feature_id = TRUE,
    unique = TRUE
  )

  expect_equal(sum(taxa$taxon_name == "Gadus" & taxa$taxon_rank == "Genus"), 1)
  expect_equal(anyDuplicated(paste(taxa$taxon_name, taxa$taxon_rank, sep = "\r")), 0)
})

test_that("extract_taxa_for_evidence works with ps_test_data_euk", {
  ps <- load_ps_test_data_euk()

  taxa <- extract_taxa_for_evidence(ps)

  expect_s3_class(taxa, "tbl_df")
  expect_true(nrow(taxa) > 0)
  expect_equal(colnames(taxa), c("taxon_name", "taxon_rank"))
})

test_that("ps_test_data_euk highest mode returns at most one taxon per feature", {
  ps <- load_ps_test_data_euk()

  taxa <- extract_taxa_for_evidence(
    ps,
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_equal(anyDuplicated(taxa$feature_id), 0)
  expect_lte(nrow(taxa), phyloseq::ntaxa(ps))
})

test_that("ps_test_data_euk default ranks exclude broad ranks", {
  ps <- load_ps_test_data_euk()

  taxa <- extract_taxa_for_evidence(ps)

  expect_false(any(taxa$taxon_rank %in% c("Domain", "Superkingdom", "Kingdom")))
  expect_true(all(taxa$taxon_rank %in% c("Phylum", "Class", "Order", "Family", "Genus", "Species")))
})

test_that("ps_test_data_euk can include feature_id for feature-level output", {
  ps <- load_ps_test_data_euk()

  taxa <- extract_taxa_for_evidence(
    ps,
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_true("feature_id" %in% colnames(taxa))
  expect_false(any(is.na(taxa$feature_id)))
})

test_that("ps_test_data_euk mode all returns multiple ranks where available", {
  ps <- load_ps_test_data_euk()

  taxa <- extract_taxa_for_evidence(
    ps,
    mode = "all",
    include_feature_id = TRUE,
    unique = FALSE
  )

  expect_gt(nrow(taxa), length(unique(taxa$feature_id)))
  expect_true(any(duplicated(taxa$feature_id)))
})

test_that("ps_test_data_euk extraction drops placeholder and empty taxa", {
  ps <- load_ps_test_data_euk()

  taxa <- extract_taxa_for_evidence(ps, mode = "all", unique = FALSE)

  expect_false(any(is.na(taxa$taxon_name)))
  expect_false(any(!nzchar(trimws(taxa$taxon_name))))
  expect_false(any(tolower(taxa$taxon_name) %in% c("uncultured", "unidentified", "unclassified", "unknown")))
  expect_false(any(grepl("^[a-z]__\\s*$", tolower(taxa$taxon_name))))
})
