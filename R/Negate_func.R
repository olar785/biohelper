#' Load a %ni%
#'
#' This function performs a Negate operation similar to %in%.
#' @export
#' @examples
#' str_a = c("A","B","C")
#' "A" %ni% str_a
#' "D" %ni% str_a
"%ni%" = Negate("%in%")
