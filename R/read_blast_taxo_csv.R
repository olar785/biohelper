#' Read raw BLAST taxonomy output
#'
#' @description
#' Reads the raw, unprocessed BLAST taxonomy output produced by
#' [blastn_taxo_assignment()], before LCA processing. The reader supports the
#' 12-column `blastn` and `megablast` output format written by that function as
#' either comma-delimited or tab-delimited text.
#'
#' @param path A single path to a raw BLAST taxonomy CSV or TSV file.
#'
#' @return A tibble with BLAST taxonomy columns named and numeric fields
#' converted to numeric or integer types.
#'
#' @export
#'
#' @examples
#' path <- tempfile(fileext = ".csv")
#' writeLines(
#'   paste(
#'     "asv1", 100, 99.5, "subject1", 0, "1e-20", 80,
#'     123, "Species name", "blast name", 98, 97,
#'     sep = ","
#'   ),
#'   path
#' )
#' read_blast_taxo_csv(path)
read_blast_taxo_csv <- function(path) {
  if (!is.character(path) || length(path) != 1 || is.na(path)) {
    stop("`path` must be a single file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  cols <- c(
    "qseqid", "qlen", "pident", "sseqid", "sgi", "evalue", "bitscore",
    "staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp"
  )

  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) == 0) {
    stop("`path` is empty.", call. = FALSE)
  }

  tab_count <- stringr::str_count(first_line, "\t")
  comma_count <- stringr::str_count(first_line, ",")
  delim <- if (tab_count >= 11 || tab_count > comma_count) "\t" else ","

  x <- readr::read_delim(
    path,
    delim = delim,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  if (ncol(x) != length(cols)) {
    stop(
      "Expected ", length(cols), " columns, but found ", ncol(x),
      ". Detected delimiter was: ", ifelse(delim == "\t", "tab", "comma"),
      call. = FALSE
    )
  }

  names(x) <- cols
  x$qlen <- as.integer(x$qlen)
  x$pident <- as.numeric(x$pident)
  x$evalue <- as.numeric(x$evalue)
  x$bitscore <- as.numeric(x$bitscore)
  x$qcovs <- as.numeric(x$qcovs)
  x$qcovhsp <- as.numeric(x$qcovhsp)

  x
}
