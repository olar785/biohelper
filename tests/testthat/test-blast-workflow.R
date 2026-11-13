test_that("all search methods return the corrected CSV and reuse every completed stage", {
  for (method in c("megablast", "blastn", "both")) {
    setup <- workflow_setup()
    setup$method <- method
    calls <- character()
    testthat::local_mocked_bindings(
      .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
      .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
        calls <<- c(calls, task)
        workflow_write_hits(workflow_hits(ids = "asv1"), output)
      }
    )
    messages <- capture.output(result <- do.call(blastn_taxo_assignment, setup), type = "message")
    expect_equal(calls, if (method == "both") c("megablast", "blastn") else method)
    expect_identical(result$asv, c("asv1", "asv2"))
    expect_named(result, c("asv", tolower(biohelper:::.biohelper_taxonomy_output_ranks()), biohelper:::.biohelper_blast_stats_columns()))
    expect_type(result$pident, "double")
    expect_type(result$qcovs, "double")
    expect_type(result$evalue, "double")
    expect_equal(result$species[[1]], "Homo sapiens")
    expect_true(is.na(result$species[[2]]))
    expect_false(file.exists(file.path(setup$output_path, "blastn_taxo_assingment.csv")))
    path <- file.path(setup$output_path, "blastn_taxo_assignment.csv")
    expect_match(readLines(path, n = 1), "asv,domain,kingdom")
    expect_equal(nrow(readr::read_csv(path, show_col_types = FALSE)), 2L)
    before <- calls
    messages <- capture.output(reused <- do.call(blastn_taxo_assignment, setup), type = "message")
    expect_identical(reused, result)
    expect_identical(calls, before)
    expect_true(any(grepl("SKIPPED (valid checkpoint)", messages, fixed = TRUE)))
    metadata <- readRDS(file.path(setup$output_path, "blastn_taxo_assignment_provenance.rds"))
    expect_true(metadata$complete)
    expect_true(metadata$input_provenance_verified)
    expect_equal(metadata$query_count, 2L)
  }
})

test_that("resumption is safe after interruption immediately following any stage", {
  original_stage <- biohelper:::.biohelper_stage
  for (interrupt_after in c("megablast_search", "blastn_search", "megablast_taxonomy", "blastn_taxonomy", "consensus")) {
    setup <- workflow_setup()
    searches <- character()
    completed <- character()
    interrupt <- TRUE
    testthat::local_mocked_bindings(
      .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
      .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
        searches <<- c(searches, task)
        workflow_write_hits(workflow_hits(), output)
      },
      .biohelper_stage = function(name, ...) {
        result <- original_stage(name, ...)
        completed <<- c(completed, name)
        if (interrupt && identical(name, interrupt_after)) stop("Simulated interruption")
        result
      }
    )
    expect_error(do.call(blastn_taxo_assignment, setup), "Simulated interruption")
    before <- searches
    interrupt <- FALSE
    result <- do.call(blastn_taxo_assignment, setup)
    expect_equal(nrow(result), 2L)
    expect_equal(sum(searches == "megablast"), 1L)
    expect_equal(sum(searches == "blastn"), 1L)
    expect_true(all(before %in% searches))
  }
})

test_that("partially written searches never replace a successful output", {
  setup <- workflow_setup()
  calls <- character()
  fail <- TRUE
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- c(calls, task)
      if (task == "blastn" && fail) { writeLines("partial", output); stop("BLAST failed halfway") }
      workflow_write_hits(workflow_hits(), output)
    }
  )
  expect_error(do.call(blastn_taxo_assignment, setup), "failed halfway")
  expect_true(file.exists(file.path(setup$output_path, "megablast_output.csv")))
  expect_false(file.exists(file.path(setup$output_path, "blastn_output.csv")))
  fail <- FALSE
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn", "blastn"))
})

test_that("taxonomy/merge failures cannot force another expensive BLAST search", {
  setup <- workflow_setup()
  calls <- character()
  original_prepare <- biohelper:::.biohelper_prepare_lca_db
  broken <- TRUE
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- c(calls, task); workflow_write_hits(workflow_hits(), output)
    },
    .biohelper_prepare_lca_db = function(...) {
      if (broken) stop("Database is temporarily inaccessible")
      original_prepare(...)
    }
  )
  expect_error(do.call(blastn_taxo_assignment, setup), "temporarily inaccessible")
  expect_equal(calls, c("megablast", "blastn"))
  broken <- FALSE
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn"))
  writeLines("corrupt", file.path(setup$output_path, "blastn_taxo_assignment.csv"))
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn"))
})

test_that("parameter-aware invalidation only reruns the affected stages", {
  setup <- workflow_setup()
  searches <- character()
  executed <- character()
  original_stage <- biohelper:::.biohelper_stage
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      searches <<- c(searches, task); workflow_write_hits(workflow_hits(), output)
    },
    .biohelper_stage = function(name, ...) {
      result <- original_stage(name, ...)
      if (result$status == "completed") executed <<- c(executed, name)
      result
    }
  )
  do.call(blastn_taxo_assignment, setup)
  executed <- character()
  setup$minSim <- 98
  do.call(blastn_taxo_assignment, setup)
  expect_equal(searches, c("megablast", "blastn"))
  expect_equal(executed, c("megablast_taxonomy", "blastn_taxonomy", "consensus"))
  executed <- character()
  setup$megablast_opts <- "-evalue 1e-10 -max_target_seqs 10"
  do.call(blastn_taxo_assignment, setup)
  expect_equal(searches, c("megablast", "blastn", "megablast"))
  expect_equal(executed, c("megablast_search", "megablast_taxonomy", "consensus"))
  expect_gt(length(list.files(file.path(setup$output_path, ".biohelper/archive"))), 0L)
  executed <- character()
  writeLines("bad partial CSV", file.path(setup$output_path, "blastn_output_processed.csv"))
  do.call(blastn_taxo_assignment, setup)
  expect_equal(searches, c("megablast", "blastn", "megablast"))
  expect_equal(executed, c("blastn_taxonomy", "consensus"))
})

test_that("legacy searches are structurally adopted and final-only legacy tables migrate", {
  setup <- workflow_setup()
  dir.create(setup$output_path)
  workflow_write_hits(workflow_hits(), file.path(setup$output_path, "megablast_output.csv"))
  workflow_write_hits(workflow_hits(), file.path(setup$output_path, "blastn_output.csv"))
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(...) stop("Legacy BLAST should be reused")
  )
  warnings <- character()
  result <- withCallingHandlers(do.call(blastn_taxo_assignment, setup), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  expect_equal(sum(grepl("adopting existing file", warnings)), 2L)
  provenance <- readRDS(file.path(setup$output_path, "blastn_taxo_assignment_provenance.rds"))
  expect_false(provenance$input_provenance_verified)

  legacy <- workflow_setup()
  dir.create(legacy$output_path)
  source <- file.path(legacy$output_path, "blastn_taxo_assingment.csv")
  utils::write.table(data.frame(asv = "asv1", domain = "Eukaryota", kingdom = "Metazoa"), source, row.names = FALSE)
  checksum <- tools::md5sum(source)
  expect_warning(adopted <- do.call(blastn_taxo_assignment, legacy), "Original parameters")
  expect_identical(unname(tools::md5sum(source)), unname(checksum))
  expect_equal(adopted$asv, c("asv1", "asv2"))
  expect_true(all(is.na(adopted$evalue)))
  expect_true(file.exists(file.path(legacy$output_path, "blastn_taxo_assignment.csv")))
  expect_identical(do.call(blastn_taxo_assignment, legacy), adopted)
})

test_that("resume/overwrite combinations protect existing user files", {
  setup <- workflow_setup()
  calls <- 0L
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- calls + 1L; workflow_write_hits(workflow_hits(), output)
    }
  )
  do.call(blastn_taxo_assignment, setup)
  setup$resume <- FALSE
  expect_error(do.call(blastn_taxo_assignment, setup), "Existing results")
  expect_equal(calls, 2L)
  setup$overwrite <- TRUE
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, 4L)
  setup$resume <- TRUE
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, 4L)
})

test_that("empty completed BLAST stages are reusable and retain all query IDs", {
  setup <- workflow_setup()
  calls <- 0L
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- calls + 1L; file.create(output)
    }
  )
  result <- do.call(blastn_taxo_assignment, setup)
  expect_identical(result$asv, c("asv1", "asv2"))
  expect_true(all(is.na(result$species)))
  expect_true(all(is.na(result$pident)))
  expect_identical(do.call(blastn_taxo_assignment, setup), result)
  expect_equal(calls, 2L)
})

test_that("a failed replacement retains the checkpoint for the existing completed search", {
  setup <- workflow_setup()
  calls <- character()
  failing <- FALSE
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- c(calls, task)
      if (failing) stop("Replacement failed")
      workflow_write_hits(workflow_hits(), output)
    }
  )
  do.call(blastn_taxo_assignment, setup)
  original <- setup$megablast_opts
  if (is.null(original)) original <- formals(blastn_taxo_assignment)$megablast_opts
  setup$megablast_opts <- "-evalue 1e-10"
  failing <- TRUE
  expect_error(do.call(blastn_taxo_assignment, setup), "Replacement failed")
  setup$megablast_opts <- original
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn", "megablast"))
})

test_that("workflow locks prevent concurrent replacement and recover same-host stale locks", {
  directory <- tempfile()
  dir.create(directory)
  release <- biohelper:::.biohelper_acquire_workflow_lock(directory)
  expect_error(biohelper:::.biohelper_acquire_workflow_lock(directory), "active or unverifiable")
  release()
  lock <- file.path(directory, ".biohelper/run.lock")
  expect_false(dir.exists(lock))
  dir.create(lock)
  saveRDS(list(pid = 1073741823L, host = unname(Sys.info()[["nodename"]]), token = "old"), file.path(lock, "owner.rds"))
  expect_message(release <- biohelper:::.biohelper_acquire_workflow_lock(directory), "Recovered")
  release()
  expect_false(dir.exists(lock))
})

test_that("merging failures and invalidated source inputs reuse only unaffected stages", {
  setup <- workflow_setup()
  calls <- character()
  fail <- TRUE
  original_format <- biohelper:::.biohelper_format_blastn_taxo_assignment_output
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- c(calls, task); workflow_write_hits(workflow_hits(), output)
    },
    .biohelper_format_blastn_taxo_assignment_output = function(out) {
      if (fail) stop("Consensus interrupted")
      original_format(out)
    }
  )
  expect_error(do.call(blastn_taxo_assignment, setup), "Consensus interrupted")
  fail <- FALSE
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn"))
  # Changed taxonomy identity invalidates taxonomy, never a BLAST search.
  Sys.setFileTime(setup$taxonomy_db, Sys.time() + 5)
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn"))
  # A changed query sequence invalidates both searches.
  writeLines(c(">asv1", "AAAAAAAA", ">asv2", "CCCCCCCC"), setup$queries)
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn", "megablast", "blastn"))
  # Corrupting only one raw output does not rerun the other method.
  writeLines("corrupted", file.path(setup$output_path, "blastn_output.csv"))
  do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, c("megablast", "blastn", "megablast", "blastn", "blastn"))
})

test_that("implementation fingerprints ignore source-reference environments", {
  first <- eval(parse(text = "function(x) { x + 1L }", keep.source = TRUE))
  second <- eval(parse(text = "function(x) { x + 1L }", keep.source = TRUE), envir = new.env())
  expect_identical(biohelper:::.biohelper_implementation_hash(list(first)), biohelper:::.biohelper_implementation_hash(list(second)))
})

test_that("valid processed legacy files are adopted while empty legacy raw files are rerun", {
  setup <- workflow_setup()
  calls <- 0L
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- calls + 1L; workflow_write_hits(workflow_hits(), output)
    }
  )
  do.call(blastn_taxo_assignment, setup)
  unlink(file.path(setup$output_path, ".biohelper/stages"), recursive = TRUE)
  warnings <- character()
  withCallingHandlers(do.call(blastn_taxo_assignment, setup), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  expect_equal(calls, 2L)
  expect_equal(sum(grepl("adopting existing file", warnings)), 4L)
  empty <- workflow_setup()
  dir.create(empty$output_path)
  file.create(file.path(empty$output_path, "megablast_output.csv"))
  file.create(file.path(empty$output_path, "blastn_output.csv"))
  do.call(blastn_taxo_assignment, empty)
  expect_equal(calls, 4L)
})

test_that("a well-formed but incomplete processed legacy table is recalculated from raw hits", {
  setup <- workflow_setup()
  calls <- 0L
  testthat::local_mocked_bindings(
    .biohelper_blast_executable = function(path) list(path = path, version = "mock BLAST 1"),
    .biohelper_run_assignment_search = function(executable, queries, db, task, options, nthreads, output, log) {
      calls <<- calls + 1L
      workflow_write_hits(workflow_hits(ids = c("asv1", "asv2")), output)
    }
  )
  do.call(blastn_taxo_assignment, setup)
  path <- file.path(setup$output_path, "megablast_output_processed.csv")
  partial <- readr::read_csv(path, show_col_types = FALSE)
  readr::write_csv(partial[1, ], path)
  unlink(file.path(setup$output_path, ".biohelper/stages/megablast_taxonomy.rds"))
  result <- do.call(blastn_taxo_assignment, setup)
  expect_equal(calls, 2L)
  expect_equal(nrow(readr::read_csv(path, show_col_types = FALSE)), 2L)
  expect_true(all(result$species == "Homo sapiens"))
  metadata <- readRDS(file.path(setup$output_path, ".biohelper/stages/megablast_taxonomy.rds"))
  expect_equal(metadata$validation$rows, 2L)
})
