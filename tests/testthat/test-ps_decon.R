make_ps_decon_fixture <- function(amplicon_type) {
  otu <- matrix(
    c(
      10, 4, 0,
      8, 0, 0,
      0, 3, 7
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("sample_1", "sample_2", "control_1"),
      c("asv_sample", "asv_shared", "asv_control")
    )
  )

  sample_df <- data.frame(
    sample_id = rownames(otu),
    amplicon_type = amplicon_type,
    row.names = rownames(otu),
    stringsAsFactors = FALSE
  )

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = FALSE),
    phyloseq::sample_data(sample_df)
  )
}

ps_decon_complete_output <- function(amplicon_type) {
  suppressMessages(
    suppressWarnings(
      capture.output(
        output <- ps_decon(
          make_ps_decon_fixture(amplicon_type),
          method = "complete_asv_removal"
        )
      )
    )
  )
  output
}

test_that("ps_decon treats amplicon_type case-insensitively", {
  lowercase <- ps_decon_complete_output(c("sample", "sample", "control"))
  capitalized <- ps_decon_complete_output(c("Sample", "Sample", "Control"))
  uppercase <- ps_decon_complete_output(c("SAMPLE", "SAMPLE", "CONTROL"))

  expect_equal(phyloseq::sample_names(capitalized), phyloseq::sample_names(lowercase))
  expect_equal(phyloseq::sample_names(uppercase), phyloseq::sample_names(lowercase))
  expect_equal(phyloseq::taxa_names(capitalized), phyloseq::taxa_names(lowercase))
  expect_equal(phyloseq::taxa_names(uppercase), phyloseq::taxa_names(lowercase))
  expect_equal(pstoveg_otu(capitalized), pstoveg_otu(lowercase))
  expect_equal(pstoveg_otu(uppercase), pstoveg_otu(lowercase))
})
