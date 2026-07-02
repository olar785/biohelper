blast_taxo_cols <- c(
  "qseqid", "qlen", "pident", "sseqid", "sgi", "evalue", "bitscore",
  "staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp"
)

write_blast_taxo_fixture <- function(values, delim = ",") {
  path <- tempfile(fileext = if (identical(delim, "\t")) ".tsv" else ".csv")
  writeLines(paste(values, collapse = delim), path)
  path
}

blastn_values <- c(
  "asv_blastn", "150", "99.3", "subject_1", "0", "1e-50", "220",
  "9606", "Homo sapiens", "primates", "98", "97"
)

megablast_values <- c(
  "asv_megablast", "145", "97.8", "subject_2", "0", "2e-40", "180",
  "9913", "Bos taurus", "even-toed ungulates", "95", "94"
)

test_that("read_blast_taxo_csv reads comma-delimited raw BLAST output", {
  path <- write_blast_taxo_fixture(blastn_values, ",")

  output <- read_blast_taxo_csv(path)

  expect_s3_class(output, "tbl_df")
  expect_equal(nrow(output), 1)
  expect_equal(output$qseqid, "asv_blastn")
  expect_equal(output$sseqid, "subject_1")
})

test_that("read_blast_taxo_csv reads tab-delimited raw BLAST output", {
  path <- write_blast_taxo_fixture(blastn_values, "\t")

  output <- read_blast_taxo_csv(path)

  expect_equal(output$qseqid, "asv_blastn")
  expect_equal(output$sscinames, "Homo sapiens")
})

test_that("read_blast_taxo_csv supports blastn and megablast raw outputs", {
  blastn_output <- read_blast_taxo_csv(write_blast_taxo_fixture(blastn_values, ","))
  megablast_output <- read_blast_taxo_csv(write_blast_taxo_fixture(megablast_values, "\t"))

  expect_equal(blastn_output$qseqid, "asv_blastn")
  expect_equal(megablast_output$qseqid, "asv_megablast")
  expect_equal(colnames(blastn_output), colnames(megablast_output))
})

test_that("read_blast_taxo_csv converts numeric columns", {
  output <- read_blast_taxo_csv(write_blast_taxo_fixture(blastn_values, ","))

  expect_type(output$qlen, "integer")
  expect_type(output$pident, "double")
  expect_type(output$evalue, "double")
  expect_type(output$bitscore, "double")
  expect_type(output$qcovs, "double")
  expect_type(output$qcovhsp, "double")
})

test_that("read_blast_taxo_csv errors when column count is not 12", {
  path <- write_blast_taxo_fixture(blastn_values[-length(blastn_values)], ",")

  expect_error(
    read_blast_taxo_csv(path),
    "Expected 12 columns, but found 11"
  )
})

test_that("read_blast_taxo_csv returns expected BLAST taxonomy columns", {
  output <- read_blast_taxo_csv(write_blast_taxo_fixture(blastn_values, ","))

  expect_equal(colnames(output), blast_taxo_cols)
})
