% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxo_bar_plot.R
\name{taxo_bar_plot}
\alias{taxo_bar_plot}
\title{taxo_bar_plot}
\usage{
taxo_bar_plot(
  ps_obj,
  rank1 = "Phylum",
  rank2 = "Family",
  n_rank1 = NA,
  n_rank2 = 6,
  colors = c("cyan", "palegreen", "yellow", "deeppink ", "white", "dodgerblue",
    "lightsalmon"),
  label = NA,
  alpha_num_ordering = F
)
}
\arguments{
\item{rank1}{Taxonomic rank to be associated with colours (e.g. Phylum)}

\item{rank2}{Taxonomic rank to be associated with shades of colours (e.g. Family)}

\item{n_rank2}{Number of shades to use per colour. Taxa which rank below that number will be assigned to "Other".}

\item{colors}{Colors to be used}

\item{alpha_num_ordering}{Whether to order samples alpha-numerically (default == F)}
}
\description{
This function takes in a phyloseq object with taxonomy and returns a bar plot with colours and shades per
taxa. The output is not directly customizable but can be re-used to add facets afterward.
Sample can also be ordered using the s_order argument.
}
\examples{
colors = c("cyan", "palegreen", "yellow", "deeppink ", "white", "dodgerblue", "lightsalmon")\cr
ps_test_data_t = ps_test_data \%>\% tax_glom('Family') \%>\% microbiome::transform(transform = "compositional")\cr
p1 = taxo_bar_plot(ps_test_data_t, rank1 = "Phylum", rank2 = "Family", colors = colors)\cr
p1 + facet_wrap(extraction_method~., drop = TRUE, scale="free", nrow = 1) + ggtitle("Taxonomic composition per extraction method")
}
