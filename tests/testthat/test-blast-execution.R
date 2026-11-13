test_that("raw BLAST reader rejects malformed rows and non-numeric statistics", {
  file <- workflow_write_hits(workflow_hits())
  lines <- readLines(file)
  writeLines(sub("99", "not_numeric", lines, fixed = TRUE), file)
  expect_error(read_blast_taxo_csv(file), "Invalid BLAST numeric")
  file <- workflow_write_hits(workflow_hits())
  writeLines(c(readLines(file), "truncated\trow"), file)
  expect_error(read_blast_taxo_csv(file), "Malformed BLAST")
  processed <- tempfile()
  writeLines(c("ASVs,taxonomy,pident,qcovs,evalue,bitscore,sseqid", "asv1"), processed)
  expect_error(biohelper:::.biohelper_read_assignment_csv(processed, processed = TRUE), "Malformed assignment CSV")
})

test_that("BLAST execution quotes paths, checks exit codes and preserves stderr", {
  testthat::skip_on_os("windows")
  directory <- tempfile("command paths ")
  dir.create(directory)
  executable <- file.path(directory, "fake blastn")
  script <- c("#!/bin/sh", 'out=""', 'while [ "$#" -gt 0 ]; do',
    '  if [ "$1" = "-out" ]; then shift; out="$1"; fi', '  shift', 'done',
    'printf "asv1\\t100\\t99\\tsubject\\t0\\t1e-30\\t200\\t9\\tHomo sapiens\\tprimates\\t100\\t100\\n" > "$out"',
    'echo "retained diagnostic message" >&2', "exit 0")
  writeLines(script, executable)
  Sys.chmod(executable, "0755")
  output <- file.path(directory, "output path.csv")
  log <- file.path(directory, "stderr log.txt")
  biohelper:::.biohelper_run_assignment_search(executable, file.path(directory, "query path.fa"),
    file.path(directory, "db prefix"), "megablast", "-evalue 1e-10", 1, output, log)
  expect_equal(read_blast_taxo_csv(output)$qseqid, "asv1")
  expect_match(paste(readLines(log), collapse = "\n"), "retained diagnostic message")
  writeLines(c("#!/bin/sh", 'echo "deliberate failure" >&2', "exit 7"), executable)
  expect_error(biohelper:::.biohelper_run_assignment_search(executable, "query", "db", "blastn", "", 1, output, log), "exit 7.*deliberate failure")
  expect_match(paste(readLines(log), collapse = "\n"), "deliberate failure")
  expect_error(biohelper:::.biohelper_run_assignment_search(executable, "query", "db", "blastn", "-out unsafe", 1, output, log), "function arguments")
  expect_error(biohelper:::.biohelper_run_assignment_search(executable, "query", "db", "blastn", "-evalue 1; echo unsafe", 1, output, log), "shell control")
})
