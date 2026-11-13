test_that("missing, inaccessible and corrupt taxonomy databases fail explicitly without downloads", {
  testthat::local_mocked_bindings(.biohelper_download_ncbi_archive = function(...) stop("Unexpected download"))
  expect_error(biohelper:::.biohelper_prepare_lca_db(FALSE, tempfile()), "database is missing")
  broken <- tempfile()
  writeLines("not a database", broken)
  expect_error(biohelper:::.biohelper_prepare_lca_db(FALSE, broken), "database.*(corrupt|functional|open)|Cannot open")
  incomplete <- tempfile()
  con <- DBI::dbConnect(RSQLite::SQLite(), incomplete)
  DBI::dbWriteTable(con, "unrelated", data.frame(x = 1))
  DBI::dbDisconnect(con)
  expect_error(biohelper:::.biohelper_prepare_lca_db(FALSE, incomplete), "Required NCBI tables")
  valid <- workflow_database()
  info <- biohelper:::.biohelper_prepare_lca_db(FALSE, valid)
  expect_equal(info$path, normalizePath(valid))
})

test_that("HTTPS staging builds a valid replacement without touching external databases", {
  archive <- testthat::test_path("fixtures", "ncbi-prunella", "taxdmp.zip")
  downloads <- 0L
  testthat::local_mocked_bindings(.biohelper_download_ncbi_archive = function(destination) {
    downloads <<- downloads + 1L
    file.copy(archive, destination)
  })
  target <- tempfile(fileext = ".sqlite")
  result <- biohelper:::.biohelper_update_ncbi_safe(target)
  expect_equal(downloads, 1L)
  expect_true(file.exists(target))
  expect_equal(result$provenance$format, "biohelper_ncbi_compact_v1")
  expect_equal(result$provenance$source, "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
  lookup <- biohelper:::.biohelper_lca_lookup(c("39176", "175132"), result)
  expect_equal(lookup$resolved, 2L)
  checksum <- unname(tools::md5sum(target))
  external <- workflow_database("taxizedb")
  external_checksum <- unname(tools::md5sum(external))
  expect_error(biohelper:::.biohelper_update_ncbi_safe(external), "externally managed")
  expect_identical(unname(tools::md5sum(external)), external_checksum)
  expect_equal(downloads, 1L)

  # A bad transfer cannot delete or alter the current working database.
  testthat::local_mocked_bindings(.biohelper_download_ncbi_archive = function(destination) writeLines("<html>FTP proxy error</html>", destination))
  expect_error(biohelper:::.biohelper_update_ncbi_safe(target), "not a valid ZIP")
  expect_identical(unname(tools::md5sum(target)), checksum)
  expect_equal(biohelper:::.biohelper_validate_ncbi_db(target)$path, normalizePath(target))
})

test_that("successful updates preserve an archived previous database and failed builds preserve it in place", {
  archive <- testthat::test_path("fixtures", "ncbi-prunella", "taxdmp.zip")
  directory <- tempfile()
  dir.create(directory)
  target <- file.path(directory, "ncbi.sqlite")
  testthat::local_mocked_bindings(.biohelper_download_ncbi_archive = function(destination) file.copy(archive, destination))
  biohelper:::.biohelper_update_ncbi_safe(target)
  checksum <- unname(tools::md5sum(target))
  biohelper:::.biohelper_update_ncbi_safe(target)
  archived <- list.files(file.path(directory, "archive"), full.names = TRUE)
  expect_equal(length(archived), 1L)
  expect_identical(unname(tools::md5sum(archived)), checksum)
  current <- unname(tools::md5sum(target))
  testthat::local_mocked_bindings(.biohelper_build_ncbi_archive = function(archive, output, workdir) {
    writeLines("interrupted build", output)
    stop("Build failed")
  })
  expect_error(biohelper:::.biohelper_update_ncbi_safe(target), "Build failed")
  expect_identical(unname(tools::md5sum(target)), current)
})

test_that("both-method workflow updates once and internal LCA calls never update", {
  setup <- workflow_setup()
  updates <- 0L
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) workflow_write_hits(workflow_hits(), output),
    .biohelper_update_ncbi_safe = function(path) {
      updates <<- updates + 1L
      biohelper:::.biohelper_validate_ncbi_db(setup$taxonomy_db)
    }
  )
  do.call(blastn_taxo_assignment, c(setup, list(update = TRUE)))
  expect_equal(updates, 1L)
  do.call(blastn_taxo_assignment, c(setup, list(update = FALSE)))
  expect_equal(updates, 1L)
  lcaPident(workflow_write_hits(workflow_hits()), tempfile(), taxonomy_db = setup$taxonomy_db, update = TRUE)
  expect_equal(updates, 2L)
})

test_that("database read permissions and cyclic ancestry are reported explicitly", {
  database <- workflow_database()
  Sys.chmod(database, "0000")
  on.exit(Sys.chmod(database, "0600"), add = TRUE)
  if (file.access(database, 4L) != 0L) {
    expect_error(biohelper:::.biohelper_validate_ncbi_db(database), "inaccessible")
  }
  Sys.chmod(database, "0600")
  connection <- DBI::dbConnect(RSQLite::SQLite(), database)
  DBI::dbExecute(connection, "UPDATE nodes SET parent=9 WHERE id=8")
  DBI::dbDisconnect(connection)
  info <- biohelper:::.biohelper_validate_ncbi_db(database)
  expect_error(biohelper:::.biohelper_lca_lookup("9", info), "cyclic or excessive")
})
