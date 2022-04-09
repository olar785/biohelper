#' Loads tidy_psmelt
#'
#'The tidy_psmelt function is a specialized melt function for melting phyloseq objects (instances of the phyloseq class), usually for producing graphics with ggplot2.
#'This function was originally created and posted by mikemc (https://github.com/mikemc) to solve an issue with phyloseq psmelt function (https://github.com/joey711/phyloseq/issues/1035)
#'
#' @export
#' @examples
#'
#' tidy_psmelt(ps_test_data)

tidy_psmelt <- function(physeq) {
  # First, make sure taxa are rows
  if (taxa_are_rows(physeq)) {
    st <- otu_table(physeq) %>% t
  } else {
    st <- otu_table(physeq)
  }
  # Convert to a tibble in tidy form
  st <- st %>%
    as("matrix") %>%
    as_tibble(rownames = "Sample") %>%
    gather("OTU", "Abundance", -Sample)
  # Get the sample data as a tibble
  sam <- sample_data(physeq) %>%
    as("data.frame") %>%
    as_tibble(rownames = "Sample")
  # Same for the tax table
  tax <- tax_table(physeq) %>%
    as("matrix") %>%
    as_tibble(rownames = "OTU")
  # Now join them all into one
  tb <- st %>%
    left_join(sam, by = "Sample") %>%
    left_join(tax, by = "OTU")
  tb
}
