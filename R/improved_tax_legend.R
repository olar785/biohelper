#' improved_tax_legend
#'
#' @description
#' This function takes in a plot from the taxo_bar_plot function of biohelper and returns the figure with a modified legend, where the higher rank (e.g. Phylum [rank1]) is displayed in bold only once
#' above its associated lower groups (e.g. Families [rank2]).
#' The function is simply meant to provide a cleaner legend for this type of figure.
#'
#' @param
#' p              Plot from the taxo_bar_plot function of biohelper
#' @param
#' rank1          Higher taxonomic rank associated with colours (e.g. Phylum)
#' @param
#' rank2          Lower taxonomic rank associated with shades of colours (e.g. Family)
#' @param
#' others_last    Optionally places the groups "Others" last, within each higher taxonomic group
#' @param
#' legend_position Position of the legend
#'
#' @export
#' @examples
#' ps_test_data_t = ps_test_data %>% tax_glom('Family') %>% subset_samples(extraction_method!="na") %>% microbiome::transform(transform = "compositional")\cr
#' p1 = taxo_bar_plot(ps_test_data_t, rank1 = "Phylum", rank2 = "Family") + facet_wrap(extraction_method~., drop = TRUE, scale="free", nrow = 1) + ggtitle("Taxonomic composition per extraction method")\cr
#' improved_tax_legend(p1,rank1 = "phylum", rank2 = "family")
#'

improved_tax_legend <- function(p, rank1 = "Phylum", rank2 = "Family", others_last = FALSE, legend_position = "bottom") {

  p = p+theme(legend.position = legend_position)
  gb <- ggplot2::ggplot_build(p)

  df <- gb$plot$data
  if (is.null(df) || !length(df)) df <- gb$data[[1]]

  # rank2 column holds the combined key "rank1.rank2" (case-insensitive)
  r2_col <- names(df)[match(tolower(rank2), tolower(names(df)))]
  if (is.na(r2_col)) stop(sprintf("Could not find a '%s' column in the plot data/layers.", rank2))

  pairs <- df |>
    dplyr::mutate(
      R1 = sub("\\..*$", "", .data[[r2_col]]),
      R2 = sub("^[^.]*\\.", "", .data[[r2_col]])
    ) |>
    dplyr::distinct(R1, R2, !!rlang::sym(r2_col))

  pairs_sorted <- pairs |>
    dplyr::mutate(is_others = if (others_last) tolower(R2) %in% c("others","other") else FALSE) |>
    dplyr::arrange(R1, is_others, tolower(R2))

  header_keys  <- paste0("HDR;;", unique(pairs_sorted$R1))
  legend_order <- unlist(lapply(split(pairs_sorted, pairs_sorted$R1), function(d) {
    c(paste0("HDR;;", d$R1[1]), d[[r2_col]])
  }))

  hdr_labels   <- stats::setNames(paste0("<b>", unique(pairs_sorted$R1), "</b>"), header_keys)
  child_labels <- stats::setNames(paste0("&nbsp;&nbsp;", pairs_sorted$R2), pairs_sorted[[r2_col]])
  legend_labels <- c(hdr_labels, child_labels)

  fill_scale <- gb$plot$scales$get_scales("fill")
  r2_levels  <- unique(pairs_sorted[[r2_col]])
  r2_cols    <- stats::setNames(fill_scale$map(r2_levels), r2_levels)

  fill_values <- c(stats::setNames(rep("#00000000", length(header_keys)), header_keys), r2_cols)

  # -------- FIX: build hdr_df with one row per header key --------
  hdr_df <- data.frame(
    x = rep(0, length(header_keys)),
    y = rep(0, length(header_keys))
  )
  hdr_df[[r2_col]] <- header_keys
  # ---------------------------------------------------------------

  p +
    ggplot2::geom_point(
      data = hdr_df,
      ggplot2::aes(x = x, y = y, fill = !!rlang::sym(r2_col)),
      inherit.aes = FALSE, alpha = 0, size = 0.01, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = legend_order,
      labels = legend_labels,
      guide = ggplot2::guide_legend(
        title = paste0(rank1, " & ", rank2),
        byrow = TRUE, label.hjust = 0
      )
    ) +
    ggplot2::theme(
      legend.text  = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.key.width  = grid::unit(0.9, "lines"),
      legend.key.height = grid::unit(0.5, "lines")
    ) +
    guides(fill = guide_legend(title = paste0( paste(toupper(substr(rank1, 1, 1)), substr(rank1, 2, nchar(rank1)), sep="")," & ",paste(toupper(substr(rank2, 1, 1)), substr(rank2, 2, nchar(rank2)), sep=""))))


}
