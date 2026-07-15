upset_to_venn_parse_description_field <- function(field) {
  if (is.na(field) || !nzchar(field)) {
    return(character())
  }

  trimws(unlist(strsplit(field, ",", fixed = TRUE), use.names = FALSE))
}

test_that("eulerr is optional package metadata", {
  desc <- read.dcf(system.file("DESCRIPTION", package = "biohelper"))[1, ]
  imports <- upset_to_venn_parse_description_field(desc[["Imports"]])
  suggests <- upset_to_venn_parse_description_field(desc[["Suggests"]])

  testthat::expect_false("eulerr" %in% imports)
  testthat::expect_true("eulerr" %in% suggests)
})

test_that("NAMESPACE does not import eulerr", {
  namespace <- readLines(system.file("NAMESPACE", package = "biohelper"))

  testthat::expect_false(any(grepl("^import\\(eulerr\\)", namespace)))
  testthat::expect_false(any(grepl("^importFrom\\(eulerr,", namespace)))
})

test_that("upset_to_venn reports that eulerr is required for this function", {
  testthat::with_mocked_bindings(
    testthat::expect_error(
      upset_to_venn(data.frame(a = c(1, 0), b = c(0, 1))),
      regexp = "eulerr.*upset_to_venn"
    ),
    check_installed = function(pkg, reason = NULL, ...) {
      testthat::expect_identical(pkg, "eulerr")
      testthat::expect_identical(
        reason,
        "to use `biohelper::upset_to_venn()`"
      )
      rlang::abort(
        paste("eulerr is required", reason)
      )
    },
    .package = "rlang"
  )
})
