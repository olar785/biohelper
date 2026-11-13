consensus_table <- function(id = "asv1", genus = "Homo", species = NA_character_, pident = 99,
                            qcovs = 100, evalue = 1e-30, bitscore = 200, sseqid = "reference") {
  data.frame(ASVs = id, domain = "Eukaryota", kingdom = "Metazoa", phylum = "Chordata",
    class = "Mammalia", order = "Primates", family = "Hominidae", genus = genus, species = species,
    pident = pident, qcovs = qcovs, evalue = evalue, bitscore = bitscore, sseqid = sseqid,
    stringsAsFactors = FALSE)
}

test_that("more resolved non-conflicting assignments keep their own alignment statistics", {
  mega <- consensus_table(species = "Homo sapiens", pident = 97, qcovs = 95, evalue = 1e-20)
  blast <- consensus_table(pident = 99, qcovs = 100, evalue = 1e-40)
  result <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = mega, blastn = blast))
  expect_equal(result$Species, "Homo sapiens")
  expect_equal(result$assignment_method, "megablast")
  expect_equal(result$support_method, "megablast")
  expect_equal(result$pident, 97)
  expect_equal(result$qcovs, 95)
  expect_equal(result$evalue, 1e-20)
  swapped <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = blast, blastn = mega))
  expect_equal(swapped$assignment_method, "blastn")
  expect_equal(swapped$support_method, "blastn")
  expect_equal(swapped$pident, 97)
})

test_that("identical classifications choose support deterministically without claiming conflict resolution", {
  mega <- consensus_table(pident = 99, evalue = 0, bitscore = 200, sseqid = "b")
  blast <- consensus_table(pident = 98, evalue = 0, bitscore = 300, sseqid = "a")
  result <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = mega, blastn = blast))
  expect_equal(result$assignment_method, "blastn")
  expect_equal(result$support_method, "blastn")
  expect_equal(result$pident, 98)
  expect_false(result$assignment_method == "consensus")
  tied <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(blastn = mega, megablast = mega))
  expect_equal(tied$assignment_method, "megablast")
  expect_equal(tied$support_method, "megablast")
})

test_that("conflicting assignments truncate at their first disagreement and keep one support hit", {
  mega <- consensus_table(genus = "Homo", species = "Homo sapiens", pident = 98, qcovs = 92, evalue = 1e-20)
  blast <- consensus_table(genus = "Pan", species = "Pan troglodytes", pident = 99, qcovs = 100, evalue = 1e-30)
  result <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = mega, blastn = blast))
  expect_equal(result$Family, "Hominidae")
  expect_true(is.na(result$Genus))
  expect_true(is.na(result$Species))
  expect_equal(result$assignment_method, "consensus")
  expect_equal(result$support_method, "blastn")
  expect_equal(unname(unlist(result[, c("pident", "qcovs", "evalue")])), c(99, 100, 1e-30))
  mega$kingdom <- "Viridiplantae"
  result <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = mega, blastn = blast))
  expect_equal(result$Domain, "Eukaryota")
  expect_true(all(is.na(result[, biohelper:::.biohelper_taxonomy_output_ranks()[-1]])))
})

test_that("single-method ASVs, missing statistics and duplicates are handled predictably", {
  mega <- rbind(consensus_table("003"), consensus_table("001"))
  blast <- rbind(consensus_table("002", pident = NA_real_, qcovs = NA_real_, evalue = NA_real_), consensus_table("003"))
  first <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = mega, blastn = blast))
  second <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(blastn = blast[2:1, ], megablast = mega[2:1, ]))
  rownames(first) <- rownames(second) <- NULL
  expect_identical(first, second)
  expect_equal(first$ASV, c("001", "002", "003"))
  expect_equal(first$assignment_method[1:2], c("megablast", "blastn"))
  expect_true(is.na(first$support_method[[2]]))
  duplicated <- biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = rbind(mega, mega)))
  expect_equal(nrow(duplicated), 2L)
  conflict <- rbind(mega, consensus_table("001", genus = "Pan"))
  expect_error(biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = conflict)), "Conflicting duplicate")
})

test_that("statistics follow species in genuine CSV without becoming taxonomy", {
  input <- consensus_table(species = 'Homo sapiens, isolate "A"')
  out <- biohelper:::.biohelper_format_blastn_taxo_assignment_output(
    biohelper:::.biohelper_merge_blast_taxonomy_tables(list(megablast = input)))
  expect_named(out, c("asv", tolower(biohelper:::.biohelper_taxonomy_output_ranks()), biohelper:::.biohelper_blast_stats_columns()))
  path <- tempfile(fileext = ".csv")
  readr::write_csv(out, path, na = "NA")
  roundtrip <- readr::read_csv(path, show_col_types = FALSE)
  expect_equal(roundtrip$species, input$species)
  expect_type(roundtrip$pident, "double")
  expect_type(roundtrip$qcovs, "double")
  expect_type(roundtrip$evalue, "double")
  standard <- biohelper:::.biohelper_standardise_taxonomy_columns(out)
  expect_false("assignment_method" %in% names(standard))
  expect_equal(standard$Domain, "Eukaryota")
})

test_that("normalisation preserves numeric/provenance metadata without looking it up as taxonomy", {
  database <- ncbi_prunella_database()
  input <- data.frame(feature_id = "asv1", a = "Viridiplantae", b = "Prunella",
    pident = 99, qcovs = 100, evalue = 1e-30, assignment_method = "megablast", support_method = "Passeroidea")
  result <- suppressWarnings(taxo_normalisation(input, database, addExtra = FALSE))
  expect_equal(result$family, "Lamiaceae")
  expect_equal(result$assignment_method, "megablast")
  expect_equal(result$support_method, "Passeroidea")
  expect_type(result$pident, "double")
  expect_equal(result$evalue, 1e-30)
})
