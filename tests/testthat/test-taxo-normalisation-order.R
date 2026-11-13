test_that("real NCBI resolution accepts any number of ordered unnamed ranks", {
  database <- ncbi_prunella_database()
  normalise <- function(x, ...) suppressWarnings(taxo_normalisation(
    x, database, addExtra = FALSE, ...
  ))
  paths <- list(
    c("Lamiaceae"),
    c("Viridiplantae", "Prunella"),
    c("Eukaryota", "Viridiplantae", "Pentapetalae", "Prunella"),
    c("Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria",
      "Deuterostomia", "Chordata", "Craniata", "Vertebrata", "Gnathostomata",
      "Tetrapoda", "Amniota", "Aves", "Passeriformes", "Passeroidea", "Prunella")
  )
  for (i in seq_along(paths)) {
    taxonomy <- matrix(paths[[i]], nrow = 1, dimnames = list("asv1", paste0("level", seq_along(paths[[i]]))))
    before <- taxonomy
    expected_family <- if (i == 4L) "Prunellidae" else "Lamiaceae"
    output <- normalise(taxonomy)
    expect_equal(output$family, expected_family)
    expect_identical(taxonomy, before)

    # Deliberately misleading labels cannot impose a rank or search boundary.
    colnames(taxonomy) <- c("kingdom", "phylum", "species", "genus", paste0("other", seq_len(ncol(taxonomy))))[seq_len(ncol(taxonomy))]
    expect_identical(normalise(taxonomy), output)
    colnames(taxonomy) <- rep("clade", ncol(taxonomy))
    expect_identical(normalise(taxonomy), output)
    colnames(taxonomy) <- rep(NA_character_, ncol(taxonomy))
    expect_identical(normalise(taxonomy), output)
  }

  # Missing columns and unrecognised finer names do not set a depth limit.
  input <- data.frame(feature_id = c("plant", "bird"),
    level1 = c("Viridiplantae", "Metazoa"), level2 = NA_character_,
    level3 = "UnknownIntermediate", level4 = "Prunella",
    level5 = "UnknownFinerName")
  expect_equal(normalise(input)$family, c("Lamiaceae", "Prunellidae"))
  names(input)[2:6] <- c("kingdom", "phylum", "taxon", "taxonomy", "species")
  expect_equal(normalise(input)$family, c("Lamiaceae", "Prunellidae"))

  for (column in c("taxon", "taxonomy")) {
    packed <- data.frame(feature_id = c("plant", "bird"),
      value = c("Eukaryota;Viridiplantae;Prunella", "Eukaryota;Metazoa;Chordata;Passeroidea;Prunella"))
    names(packed)[2] <- column
    expect_equal(normalise(packed)$family, c("Lamiaceae", "Prunellidae"))
  }
})

test_that("ancestry order distinguishes candidates containing the same names", {
  lineages <- list(
    `11` = c(genus = "SharedName", clade = "GroupA", clade.1 = "GroupB", domain = "Eukaryota"),
    `22` = c(genus = "SharedName", clade = "GroupB", clade.1 = "GroupA", domain = "Eukaryota")
  )
  reference <- data.frame(family = c("Family11", "Family22"), row.names = c("11", "22"))
  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) rep("11,22", length(taxa)),
    getRawTaxonomy = function(ids, sqlFile) lineages[ids],
    getTaxonomy = function(ids, sqlFile, desiredTaxa) reference[ids, desiredTaxa, drop = FALSE],
    .package = "taxonomizr"
  )
  input <- data.frame(feature_id = c("A_then_B", "B_then_A"),
    level1 = "Eukaryota", level2 = c("GroupA", "GroupB"),
    level3 = c("GroupB", "GroupA"), level4 = "SharedName")
  output <- taxo_normalisation(input, "mock.sql", addExtra = FALSE, ranks = "family")
  expect_equal(output$family, c("Family22", "Family11"))

  # A repeated name may represent the same NCBI ancestor more than once.
  input <- data.frame(feature_id = "repeated", a = "Eukaryota", b = "GroupA",
                      c = "GroupB", d = "GroupB", e = "SharedName")
  output <- taxo_normalisation(input, "mock.sql", addExtra = FALSE, ranks = "family")
  expect_equal(output$family, "Family22")
})

test_that("column names cannot alter name cleaning or explicit epithet joining", {
  reference <- data.frame(species = "Long taxon name", row.names = "11")
  seen <- character()
  testthat::local_mocked_bindings(
    getId = function(taxa, sqlFile, onlyScientific = TRUE) {
      seen <<- c(seen, taxa)
      rep("11", length(taxa))
    },
    getTaxonomy = function(ids, sqlFile, desiredTaxa) reference[ids, desiredTaxa, drop = FALSE],
    .package = "taxonomizr"
  )
  input <- data.frame(feature_id = "asv1", kingdom = "Parent", species = "Long_taxon_name")
  output <- taxo_normalisation(input, "mock.sql", addExtra = FALSE, ranks = "species")
  expect_equal(seen, "Long taxon name")
  names(input)[2:3] <- c("level1", "level2")
  expect_identical(taxo_normalisation(input, "mock.sql", addExtra = FALSE, ranks = "species"), output)

  input <- data.frame(feature_id = "asv1", a = "Parent", b = "Long", c = "taxon_name")
  seen <- character()
  taxo_normalisation(input, "mock.sql", addExtra = FALSE, ranks = "species", spnc = TRUE)
  expect_equal(seen, "Long taxon name")
  expect_error(taxo_normalisation(input[, c("feature_id", "c")], "mock.sql", spnc = TRUE), "at least two")
})

test_that("ambiguity methods and cached BLAST are independent of input labels", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", a = "Eukaryota", b = "Prunella")
  renamed <- input
  names(renamed)[2:3] <- c("species", "phylum")
  hits <- data.frame(qseqid = "asv1", sseqid = c("record1", "record2"), staxids = "175132",
                     pident = 99, qcovhsp = 100, evalue = 1e-30, bitscore = 100)
  for (method in c("first", "lca", "unassigned", "blast_check")) {
    normalise <- function(x) suppressWarnings(taxo_normalisation(
      x, database, addExtra = FALSE, ambiguity = method, return_details = TRUE,
      blast_hits = hits, blast_options = list(min_taxa = 1L, fallback = "lca")
    ))
    expect_identical(normalise(input), normalise(renamed))
  }
})
