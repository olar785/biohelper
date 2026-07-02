#' Loads %ni%
#'
#' This function performs a Negate operation similar to %in%.
#'
#' @param x Values to match.
#' @param table Values to be matched against.
#'
#' @export
#' @examples
#' str_a = c("A","B","C")
#' "A" %ni% str_a
#' "D" %ni% str_a
"%ni%" = function(x, table) {
  !(x %in% table)
}
