barplot_test_phyloseq <- function() {
  taxonomy <- cbind(
    Domain = "Bacteria",
    Phylum = rep(c("P1", "P2", "P3"), each = 7),
    Family = paste0("F", seq_len(21))
  )
  rownames(taxonomy) <- paste0("otu", seq_len(nrow(taxonomy)))
  abundance <- matrix(seq_len(21) / sum(seq_len(21)), ncol = 1,
                      dimnames = list(rownames(taxonomy), "sample1"))
  phyloseq::phyloseq(
    phyloseq::otu_table(abundance, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(data.frame(group = "A", row.names = "sample1"))
  )
}

test_that("taxo_bar_plot preserves explicit n_rank1 and derives its default", {
  ps <- barplot_test_phyloseq()
  colours <- c("red", "blue", "green")
  for (n in c(1L, 2L)) {
    plot <- suppressMessages(taxo_bar_plot(
      ps, n_rank1 = n, n_rank2 = 2, colors = colours
    ))
    expect_s3_class(plot, "ggplot")
    phyla <- unique(sub("\\..*$", "", plot$data$Family))
    expect_equal(length(setdiff(phyla, "Others")), n)
    expect_equal(sum(plot$data$Abund), 100)
  }
  default_plot <- suppressMessages(taxo_bar_plot(ps, n_rank2 = 2, colors = colours))
  phyla <- unique(sub("\\..*$", "", default_plot$data$Family))
  expect_equal(length(setdiff(phyla, "Others")), 2L)
})

test_that("taxo_bar_plot applies n_rank2 within each rank1 group", {
  ps <- barplot_test_phyloseq()
  original <- ps
  for (n in c(1L, 3L, 6L)) {
    plot <- suppressMessages(taxo_bar_plot(
      ps, n_rank1 = 3, n_rank2 = n, colors = c("red", "blue", "green", "grey")
    ))
    for (phylum in c("P1", "P2", "P3")) {
      families <- plot$data$Family[startsWith(plot$data$Family, paste0(phylum, "."))]
      expect_equal(length(setdiff(families, paste0(phylum, ".Others"))), n)
      expect_true(paste0(phylum, ".Others") %in% families)
    }
    expect_equal(sum(plot$data$Abund), 100)
    expect_no_error(ggplot2::ggplot_build(plot))
  }
  expect_identical(ps, original)
})
