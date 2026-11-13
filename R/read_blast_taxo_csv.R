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

  x <- suppressWarnings(readr::read_delim(
    path,
    delim = delim,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    trim_ws = TRUE
  ))

  if (ncol(x) != length(cols)) {
    stop(
      "Expected ", length(cols), " columns, but found ", ncol(x),
      ". Detected delimiter was: ", ifelse(delim == "\t", "tab", "comma"),
      call. = FALSE
    )
  }

  names(x) <- cols
  if (nrow(readr::problems(x))) stop("Malformed BLAST records: inconsistent field counts or quoting in ", path, call. = FALSE)
  for (column in c("qlen", "pident", "evalue", "bitscore", "qcovs", "qcovhsp")) {
    values <- suppressWarnings(as.numeric(x[[column]]))
    invalid <- !is.finite(values) | values < 0
    if (column == "qlen") invalid <- invalid | values == 0 | values != floor(values)
    if (column %in% c("pident", "qcovs", "qcovhsp")) invalid <- invalid | values > 100
    if (any(invalid)) stop("Invalid BLAST numeric values in column '", column, "' in ", path, call. = FALSE)
    x[[column]] <- if (column == "qlen") as.integer(values) else values
  }
  if (anyNA(x$qseqid) || any(!nzchar(x$qseqid)) || anyNA(x$sseqid) || any(!nzchar(x$sseqid))) stop("BLAST query and subject IDs must be non-empty.", call. = FALSE)

  x
}
