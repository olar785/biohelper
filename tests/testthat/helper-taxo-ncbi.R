ncbi_prunella_database <- function() {
  database <- tempfile(fileext = ".sqlite")
  fixture <- testthat::test_path("fixtures", "ncbi-prunella")
  taxonomizr::read.names.sql(file.path(fixture, "names.dmp"), database)
  taxonomizr::read.nodes.sql(file.path(fixture, "nodes.dmp"), database)
  database
}
