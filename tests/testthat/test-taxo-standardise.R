test_that("old NCBI-style superkingdom is standardised to Domain", {
  input <- data.frame(
    ASV = "asv1",
    superkingdom = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    class = "Actinopteri",
    order = "Carangiformes",
    family = "Carangidae",
    genus = "Seriola",
    species = NA_character_,
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_standardise_taxonomy_columns(
    input,
    id_col = "ASV",
    output_id_col = "ASV"
  )

  expect_equal(output$Domain, "Eukaryota")
  expect_equal(output$Kingdom, "Metazoa")
  expect_false(any(tolower(colnames(output)) == "superkingdom"))
})

test_that("new NCBI-style domain is standardised to Domain", {
  input <- data.frame(
    ASV = "asv1",
    domain = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_standardise_taxonomy_columns(
    input,
    id_col = "ASV",
    output_id_col = "ASV"
  )

  expect_equal(output$Domain, "Eukaryota")
  expect_equal(output$Kingdom, "Metazoa")
  expect_false(any(tolower(colnames(output)) == "superkingdom"))
})

test_that("Domain values are preferred over legacy highest-rank values when both are present", {
  input <- data.frame(
    ASV = c("asv1", "asv2"),
    Domain = c("Eukaryota", NA_character_),
    superkingdom = c("Wrong fallback", "Bacteria"),
    Kingdom = c("Metazoa", NA_character_),
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_standardise_taxonomy_columns(
    input,
    id_col = "ASV",
    output_id_col = "ASV"
  )

  expect_equal(output$Domain, c("Eukaryota", "Bacteria"))
  expect_false(any(tolower(colnames(output)) == "superkingdom"))
})

test_that("BLAST method/source stays separate from taxonomy ranks while merging", {
  megablast <- data.frame(
    ASVs = "asv1",
    assignment_method = "megablast",
    domain = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    class = "Actinopteri",
    order = "Carangiformes",
    family = "Carangidae",
    genus = "Seriola",
    species = NA_character_,
    stringsAsFactors = FALSE
  )
  blastn <- megablast
  blastn$assignment_method <- "blastn"

  output <- biohelper:::.biohelper_merge_blast_taxonomy_tables(
    list(megablast = megablast, blastn = blastn)
  )

  expect_false("assignment_method" %in% colnames(output))
  expect_equal(output$Domain, "Eukaryota")
  expect_equal(output$Kingdom, "Metazoa")
  expect_false(identical(output$Domain, "megablast"))
})

test_that("blastn_taxo_assignment final formatter keeps method out of Domain", {
  megablast <- data.frame(
    ASVs = "0009",
    assignment_method = "megablast",
    domain = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    class = "Mammalia",
    order = "Primates",
    family = "Hominidae",
    genus = "Homo",
    species = "Homo sapiens",
    stringsAsFactors = FALSE
  )
  blastn <- megablast
  blastn$assignment_method <- "blastn"

  merged <- biohelper:::.biohelper_merge_blast_taxonomy_tables(
    list(megablast = megablast, blastn = blastn)
  )
  output <- biohelper:::.biohelper_format_blastn_taxo_assignment_output(merged)

  expect_equal(
    colnames(output),
    c("asv", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  )
  expect_false("assignment_method" %in% colnames(output))
  expect_false("superkingdom" %in% colnames(output))
  expect_equal(output$asv, "0009")
  expect_equal(output$domain, "Eukaryota")
  expect_equal(output$kingdom, "Metazoa")
  expect_equal(output$phylum, "Chordata")
  expect_equal(output$class, "Mammalia")
  expect_equal(output$order, "Primates")
  expect_equal(output$family, "Hominidae")
  expect_equal(output$genus, "Homo")
  expect_equal(output$species, "Homo sapiens")
  expect_false(any(output$domain %in% c("blastn", "megablast"), na.rm = TRUE))
})

test_that("blastn_taxo_assignment final formatter accepts legacy superkingdom input", {
  input <- data.frame(
    ASVs = "0009",
    assignment_method = "megablast",
    superkingdom = "Eukaryota",
    kingdom = "Metazoa",
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_format_blastn_taxo_assignment_output(input)

  expect_equal(output$asv, "0009")
  expect_equal(output$domain, "Eukaryota")
  expect_equal(output$kingdom, "Metazoa")
  expect_false(any(c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %in% colnames(output)))
  expect_false("superkingdom" %in% colnames(output))
  expect_false(any(output$domain %in% c("blastn", "megablast"), na.rm = TRUE))
})

test_that("blastn_taxo_assignment final formatter renames feature_id to asv", {
  input <- data.frame(
    feature_id = "0009",
    domain = "Eukaryota",
    kingdom = "Metazoa",
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_format_blastn_taxo_assignment_output(input)

  expect_equal(colnames(output)[[1]], "asv")
  expect_equal(output$asv, "0009")
  expect_equal(output$domain, "Eukaryota")
})

test_that("blastn_taxo_assignment refuses to write shifted method values in Domain", {
  shifted <- data.frame(
    ASV = "0009",
    Domain = "megablast",
    superkingdom = "Eukaryota",
    Kingdom = NA_character_,
    Phylum = "Metazoa",
    stringsAsFactors = FALSE
  )

  expect_error(
    biohelper:::.biohelper_format_blastn_taxo_assignment_output(shifted),
    "Domain contains assignment method values"
  )
  expect_error(
    biohelper:::.biohelper_validate_blastn_taxo_assignment_output(shifted),
    "Domain contains assignment method values"
  )
})

test_that("missing taxonomy keeps identifiers and missing rank columns", {
  input <- data.frame(
    ASVs = "asv1",
    taxonomy = "Unknown",
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_standardise_taxonomy_columns(
    input,
    id_col = "ASVs",
    output_id_col = "ASV"
  )

  expect_equal(output$ASV, "asv1")
  expect_true(all(is.na(output[, biohelper:::.biohelper_taxonomy_output_ranks()])))
})

test_that("regression: megablast is never shifted into Domain", {
  input <- data.frame(
    ASVs = "001decdcb27a0cb0921ad6968214737e",
    assignment_method = "megablast",
    domain = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    class = "Actinopteri",
    order = "Carangiformes",
    family = "Carangidae",
    genus = "Seriola",
    species = NA_character_,
    stringsAsFactors = FALSE
  )

  output <- biohelper:::.biohelper_merge_blast_taxonomy_tables(
    list(megablast = input)
  )

  expect_equal(output$ASV, "001decdcb27a0cb0921ad6968214737e")
  expect_equal(output$Domain, "Eukaryota")
  expect_equal(output$Kingdom, "Metazoa")
  expect_equal(output$Phylum, "Chordata")
  expect_equal(output$Class, "Actinopteri")
  expect_equal(output$Order, "Carangiformes")
  expect_equal(output$Family, "Carangidae")
  expect_equal(output$Genus, "Seriola")
  expect_false(identical(output$Domain, "megablast"))
})

test_that("taxo_merge standardises Domain without live NCBI lookup", {
  input <- data.frame(
    ASV = "asv1",
    domain = "Eukaryota",
    kingdom = "Metazoa",
    phylum = "Chordata",
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    taxo_normalisation = function(obj, sqlFile, ranks, addExtra, spnc) {
      expect_true("domain" %in% ranks)
      expect_false("superkingdom" %in% ranks)
      obj
    }
  )

  output <- taxo_merge(
    df_list = list(input),
    sqlFile = "mock.sql",
    ranks = c("Domain", "Kingdom", "Phylum")
  )

  expect_equal(output$feature_id, "asv1")
  expect_equal(output$Domain, "Eukaryota")
  expect_equal(output$Kingdom, "Metazoa")
  expect_equal(output$Phylum, "Chordata")
  expect_false(any(tolower(colnames(output)) == "superkingdom"))
})
