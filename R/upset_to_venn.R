#' upset_to_venn
#'
#' @description
#' This function takes in the dataframe created from the psto_upset function and outputs a proportional venn diagram using the
#' eulerr R package.
#' If you use this function for a publication, please cite the creators of the UpSetR and eulerr R packages (citation('UpSetR') and citation('eulerr')).
#'
#' @param
#' upset_df       Dataframe from the list returned by psto_upset. The dataframe should contain groups as columns and variables as row, using a presence (1) / absence (0) format.
#' @param
#' set_labels     Labels for groups
#' @param
#' palette        Palette of colors to use
#' @param
#' alpha          Transparency level
#' @param
#' legend         Whether to include a legend or not
#' @param
#' quant_cex      Quantity text size
#' @param
#' set_cex        Size of set names
#' @param
#' two_line       Whether to display quantity and proportion on the same line or not
#' @param
#' min_pct_to_label Hide labels below this % of union (e.g., 1.0)
#' @param
#' percent_digits Number of percent digits to use
#' @param
#' bold_sets      Set labels in bold
#'
#' @export
#' @examples
#' upset_list = psto_upset(pst = ps_test_data, grp = "extraction_method")
#' upset_to_venn(upset_df = upset_list$upset_df)
#'
upset_to_venn<- function(upset_df,
                       set_labels = NULL,
                       palette = NULL,
                       alpha = 0.75,
                       legend = TRUE,
                       # overlap controls:
                       quant_cex = 0.75,          # shrink quantity text
                       set_cex   = 1.0,           # size of set names
                       two_line  = FALSE,         # "123\n(4.7%)" instead of "123 (4.7%)"
                       min_pct_to_label = 0,      # hide labels below this % of union (e.g., 1.0)
                       percent_digits = 1,
                       bold_sets = TRUE) {

  vd = eulerr::euler(upset_df)

  # colours
  n_sets <- if (!is.null(vd$ellipses)) nrow(vd$ellipses) else 3L
  fills  <- if (is.null(palette)) viridis::viridis_pal(option = "viridis")(n_sets) else palette

  # set names styling
  lab_arg <- if (is.null(set_labels)) {
    list(cex = set_cex, font = if (bold_sets) 2 else 1)
  } else {
    list(labels = set_labels, cex = set_cex, font = if (bold_sets) 2 else 1)
  }

  if (utils::packageVersion("eulerr") >= "7.0.0" && min_pct_to_label <= 0) {
    # New API: let eulerr print counts + percents, just shrink text
    graphics::plot(
      vd,
      labels     = lab_arg,
      quantities = list(
        type = c("counts", "percent"),   # ← counts and %
        cex  = quant_cex,
        digits = percent_digits
      ),
      fill   = fills,
      alpha  = alpha,
      legend = legend
    )
  } else {
    # Fallback OR when you want to drop tiny labels via min_pct_to_label
    cnt <- vd$original.values
    pct <- 100 * cnt / sum(cnt)
    qlab <- if (two_line)
      paste0(scales::comma(cnt), "\n(", round(pct, percent_digits), "%)")
    else
      paste0(scales::comma(cnt), " (", round(pct, percent_digits), "%)")

    if (min_pct_to_label > 0) qlab[pct < min_pct_to_label] <- ""

    # When passing a character vector for quantities, we can't set cex there,
    # so set a temporary global cex for the plot call only:
    op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
    graphics::par(cex = quant_cex)

    graphics::plot(
      vd,
      labels     = lab_arg,
      quantities = qlab,   # pre-formatted "count (p%)" or blank for tiny regions
      fill   = fills,
      alpha  = alpha,
      legend = legend
    )
  }

}
