#' Resumable BLAST taxonomy assignment with provenance
#'
#' @description
#' Run local megablast, blastn or both, apply [lcaPident()], and combine compatible
#' assignments. Each search, each taxonomy table and the final consensus has an
#' independent, validated checkpoint. A failed downstream stage does not repeat
#' a completed BLAST search. Existing taxonomic thresholds and search-option
#' defaults are retained.
#'
#' @param blastapp_path Path to the BLAST+ `blastn` executable, or a command
#' resolvable on PATH. The executable and version are recorded.
#' @param method `"megablast"` (default), `"blastn"` or `"both"`.
#' @param queries Query FASTA file with unique, non-empty feature IDs. Headers
#' are read up to their first whitespace. Every query is retained in the final
#' table, including queries without qualifying hits.
#' @param megablast_opts Additional BLAST option string. Default
#' `"-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8"`, retained for compatibility.
#' BLAST interprets `-perc_identity` in percent, not as a fraction. Controlled
#' input/output/task/thread flags and shell control characters are rejected.
#' @param blastn_opts Additional blastn options. Default
#' `"-evalue 0.001 -max_target_seqs 5 -perc_identity 0.7"`, retained for compatibility.
#' @param db Formatted local BLAST database prefix or name resolvable through
#' BLASTDB. Reference identity records available database files, sizes and
#' modification times, avoiding full hashing of very large sequence databases.
#' @param output_path Output directory. Names and paths may contain spaces.
#' @param nthreads Positive integer number of BLAST threads.
#' @param minSim Species identity threshold (percent), default 97. This is not
#' an overall hit filter; see `pident` and [lcaPident()].
#' @param minCov Minimum `qcovs` query coverage (percent), default 80.
#' @param update Logical, default FALSE. TRUE builds/validates a fresh local
#' taxonomy database once, before assignments, using HTTPS and a separate
#' biohelper cache. Internal LCA calls always receive update=FALSE. No download
#' occurs when FALSE. See Database handling below.
#' @param pident `"before"` (default), `"after"` or `"no"`: apply rank identity
#' thresholds before LCA, after LCA using the representative hit, or not at all.
#' @param pgenus Genus identity threshold, default 95 percent.
#' @param pfamily Family identity threshold, default 87 percent.
#' @param porder Order identity threshold, default 83 percent.
#' @param pclass Class identity threshold, default 81 percent.
#' @param pphylum Phylum identity threshold, default 79 percent.
#' @param pkingdom Kingdom identity threshold, default 71 percent.
#' @param taxonly Retained for compatibility (including the historical string
#' `"TRUE"`). Representative statistics are now always retained in taxonomy
#' outputs; this argument does not remove them.
#' @param resume Logical, default TRUE. Reuse structurally valid, unchanged
#' completed stages. Missing, interrupted, corrupted or invalidated stages and
#' dependent stages are recalculated. Valid older files can be adopted with an
#' explicit warning that their original parameters/reference versions are unknown.
#' @param overwrite Logical, default FALSE. FALSE preserves superseded files
#' under `.biohelper/archive`, with checkpoint provenance when available.
#' TRUE allows them to be replaced after a valid replacement has completed.
#' This does not force rerunning a valid checkpoint when resume=TRUE.
#' @param taxonomy_db Optional local NCBI taxonomy SQLite file. Supports both
#' taxizedb names/nodes schemas and taxonomizr `accessionTaxa.sql` schemas.
#' Default NULL uses biohelper's managed cache if present, otherwise taxizedb's
#' existing local `NCBI.sql`. This is distinct from the BLAST sequence `db`.
#'
#' @section Checkpoints and interruption:
#' Search outputs are `megablast_output.csv` and/or `blastn_output.csv`: despite
#' these legacy suffixes, raw BLAST files remain headerless 12-column TSV.
#' Processed tables are `megablast_output_processed.csv` and/or
#' `blastn_output_processed.csv`, written as genuine CSV. The final genuine CSV
#' is `blastn_taxo_assignment.csv` (corrected spelling).
#'
#' Stage records are small RDS files under `.biohelper/stages`. They contain
#' completion status, input/output hashes, relevant parameters, software and
#' database identity, output validation and assignment summaries where available.
#' Temporary outputs are validated before publication. Atomic replacement and
#' hard-link preservation are used where supported; otherwise previous files
#' are preserved by renaming. Search stdout/stderr is retained in per-method logs.
#'
#' Rerun the same call after interruption. Changing LCA parameters invalidates
#' processed tables and consensus, not BLAST. Changing one method's search
#' options invalidates that search and its downstream results independently.
#' Corrupt final output only requires merging again. Failed replacements retain
#' the checkpoint of the completed file still present at its original path.
#' A directory lock prevents concurrent runs writing the same outputs; stale
#' same-host locks are recovered automatically when their process no longer exists.
#' Cross-host/unverifiable locks require confirming the previous process stopped.
#'
#' The resume/overwrite combinations are:
#' \describe{
#' \item{TRUE / FALSE (default)}{Reuse valid checkpoints; preserve superseded files.}
#' \item{TRUE / TRUE}{Reuse valid checkpoints; discard superseded files only on successful replacement.}
#' \item{FALSE / FALSE}{Start fresh only if requested output files/checkpoints do not already exist.}
#' \item{FALSE / TRUE}{Explicit fresh run; rerun requested stages and replace superseded results safely.}
#' }
#' Archive directories can grow across parameter changes; old versions are not
#' automatically deleted. No duplicate copy of a current large raw file is made
#' solely for checkpointing.
#'
#' @section Legacy adoption:
#' File existence alone is not sufficient. Raw legacy files require valid records
#' with expected query IDs; an empty legacy file without a checkpoint is rerun,
#' because zero hits cannot be distinguished from interruption. Completed empty
#' outputs with valid checkpoints are reusable. Legacy processed files without
#' the new representative-hit columns are reprocessed from reusable raw BLAST.
#'
#' The old filename `blastn_taxo_assingment.csv` is detected and preserved. When
#' only a structurally valid final table exists, it can be migrated without
#' rerunning BLAST. Unknown statistics/provenance remain NA and adoption is marked
#' unverified; old parameters cannot be inferred. When intermediates exist, they
#' are used to reconstruct the new statistics and the corrected final filename.
#'
#' @section Database handling:
#' Existing databases are validated read-only before taxonomy processing.
#' Missing files, permission problems, malformed schemas and SQL/ancestry failures
#' raise errors, rather than becoming silent missing assignments. Genuine missing
#' TaxIDs remain unresolved. Lookups use one read-only SQLite connection per LCA
#' call, batched unique TaxIDs and cached ancestor records; no online taxonomy
#' queries are used. Concurrent changes to a taxonomy database cause a retryable
#' failure instead of mixing database versions.
#'
#' update=TRUE downloads the official NCBI taxonomy ZIP via HTTPS into a staging
#' directory, builds compact names/nodes tables, and checks SQLite and parent/name
#' integrity before publishing. The old working cache is archived only after the
#' replacement passes. Installed taxizedb code and its shared database are not
#' modified. An explicitly supplied externally managed database is not replaced;
#' use taxonomy_db=NULL or a new path for an update. This compact cache intentionally
#' omits taxizedb's large hierarchy table and is consumed by biohelper's own reader.
#'
#' @section Consensus and representative statistics:
#' Compatible assignments retain the result with the most populated taxonomy
#' ranks. If both methods disagree, retain the selected higher-level taxonomy
#' and clear the first conflicting rank and all finer ranks. ASVs present in one
#' method are retained. Conflicting duplicate ASVs within one method are rejected.
#'
#' A representative hit is ordered by lowest E-value, highest bitscore, highest
#' identity, highest qcovs, then subject ID. Complete ties have a canonical row
#' tie-breaker. Between otherwise tied methods megablast is preferred. For equally
#' resolved compatible/identical classifications, this representative ranking
#' chooses the method; assignment_method records that method, not consensus.
#' Consensus is reserved for taxonomy actually truncated because of disagreement.
#'
#' More resolved compatible assignments retain their own method's statistics.
#' For disagreement, statistics come from the better representative method. Identity,
#' coverage and E-value always come from one hit. These are alignment measurements,
#' not LCA confidence, and need not individually justify the whole consensus.
#'
#' @return Invisibly returns the final data.frame and writes it as CSV, with
#' columns in this order: `asv`, `domain`, `kingdom`, `phylum`, `class`, `order`,
#' `family`, `genus`, `species`, `pident`, `qcovs`, `evalue`, `assignment_method`,
#' `support_method`. One row per query ASV, sorted by ID; missing values are NA
#' (written as `NA`), with no row-name column. Alignment statistics are numeric.
#' phyloseq tax_table may coerce them to character alongside rank names.
#'
#' assignment_method is megablast/blastn for a selected compatible assignment,
#' or consensus for conflict truncation. support_method identifies the one method
#' providing representative statistics; it is NA when no statistics are available.
#' ASVs without qualifying hits have missing taxonomy/statistics/provenance.
#'
#' `blastn_taxo_assignment_provenance.rds` records parameters, versions, database
#' identity, query/hit counts, stage records, unresolved TaxID counts and assignment
#' summaries. A legacy final-only adoption records its unverified provenance in
#' the consensus checkpoint. Existing callers selecting rank columns by name
#' remain compatible; callers assuming exactly nine columns must account for the
#' five appended columns. The previous implementation documented a returned table
#' but actually returned NULL; the function now returns its documented data.frame.
#'
#' @export
#' @examples
#' \dontrun{
#' arguments <- list(
#'   blastapp_path = "/tools/blastn", method = "both", queries = "ASVs.fasta",
#'   db = "/databases/nt", output_path = "blast_results", nthreads = 4,
#'   taxonomy_db = "accessionTaxa.sql"
#' )
#' result <- do.call(blastn_taxo_assignment, arguments)
#' # After interruption, use the same call: completed stages are skipped.
#' result <- do.call(blastn_taxo_assignment, arguments)
#' # Revise taxonomy thresholds without repeating the searches.
#' arguments$minSim <- 98
#' revised <- do.call(blastn_taxo_assignment, arguments)
#' # Explicit fresh run, allowing replacement after outputs validate.
#' arguments$resume <- FALSE
#' arguments$overwrite <- TRUE
#' fresh <- do.call(blastn_taxo_assignment, arguments)
#' }
#'

blastn_taxo_assignment <- function(
    blastapp_path,
    method = "megablast",
    queries,
    megablast_opts = "-evalue 0.001 -max_target_seqs 5 -perc_identity 0.8",
    blastn_opts = "-evalue 0.001 -max_target_seqs 5 -perc_identity 0.7",
    db,
    output_path,
    nthreads,
    minSim = 97,
    minCov = 80,
    update = FALSE,
    pident = "before",
    pgenus = 95,
    pfamily = 87,
    porder = 83,
    pclass = 81,
    pphylum = 79,
    pkingdom = 71,
    taxonly = "TRUE",
    resume = TRUE,
    overwrite = FALSE,
    taxonomy_db = NULL
) {
  method <- match.arg(method, c("megablast", "blastn", "both"))
  pident <- match.arg(pident, c("before", "after", "no"))
  for (flag in list(resume, overwrite, update)) if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("resume, overwrite and update must be TRUE or FALSE.", call. = FALSE)
  for (value in list(queries, db, output_path, blastapp_path)) if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) stop("Supply single paths/names for queries, db, output_path and blastapp_path.", call. = FALSE)
  if (!is.numeric(nthreads) || length(nthreads) != 1L || !is.finite(nthreads) || nthreads < 1 || nthreads != floor(nthreads)) stop("nthreads must be a positive integer.", call. = FALSE)
  nthreads <- as.integer(nthreads)
  for (value in list(minSim, minCov, pgenus, pfamily, porder, pclass, pphylum, pkingdom)) if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value > 100) stop("Identity/coverage thresholds must be between 0 and 100.", call. = FALSE)
  if (!file.exists(queries)) stop("Query FASTA file does not exist: ", queries, call. = FALSE)
  ids <- .biohelper_query_ids(queries)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_path, mustWork = TRUE)
  if (file.access(output_path, 2L) != 0L) stop("Output directory is not writable: ", output_path, call. = FALSE)
  release_lock <- .biohelper_acquire_workflow_lock(output_path)
  on.exit(release_lock())
  checkpoints <- file.path(output_path, ".biohelper", "stages")
  methods <- if (method == "both") c("megablast", "blastn") else method
  outputs <- stats::setNames(file.path(output_path, paste0(methods, "_output.csv")), methods)
  processed <- stats::setNames(file.path(output_path, paste0(methods, "_output_processed.csv")), methods)
  final_file <- file.path(output_path, "blastn_taxo_assignment.csv")
  legacy_file <- file.path(output_path, "blastn_taxo_assingment.csv")
  if (!resume && !overwrite && any(file.exists(c(outputs, processed, final_file, legacy_file)))) stop("Existing results found. Use resume=TRUE, a new output directory, or overwrite=TRUE for an explicit fresh run.", call. = FALSE)
  versions <- .biohelper_blast_versions()
  parameters <- list(minSim = minSim, minCov = minCov, pident = pident,
                     pgenus = pgenus, pfamily = pfamily, porder = porder,
                     pclass = pclass, pphylum = pphylum, pkingdom = pkingdom, taxonly = taxonly)
  options <- list(megablast = megablast_opts, blastn = blastn_opts)
  for (search in methods) .biohelper_validate_assignment_options(options[[search]])
  query_fingerprint <- .biohelper_file_fingerprint(queries)
  reference <- .biohelper_blast_db_identity(db)
  configuration <- list(method = method, queries = query_fingerprint, database = reference,
    options = options[methods], parameters = parameters, versions = versions,
    taxonomy = .biohelper_file_fingerprint(.biohelper_lca_db_path(taxonomy_db), content = FALSE))
  final_metadata_file <- file.path(checkpoints, "consensus.rds")
  final_metadata <- if (file.exists(final_metadata_file)) tryCatch(readRDS(final_metadata_file), error = function(e) NULL) else NULL
  if (resume && !update && isTRUE(final_metadata$adopted_final) &&
      identical(final_metadata$key, .biohelper_object_hash(configuration)) &&
      identical(final_metadata$output, .biohelper_file_fingerprint(final_file))) {
    out <- .biohelper_read_assignment_csv(final_file, ids = ids)
    message("Consensus merging: SKIPPED (adopted final; provenance unverified)")
    return(invisible(out))
  }
  # A final-only legacy analysis can be recovered without rerunning BLAST.
  # When intermediates exist, use them to recover the new alignment statistics.
  if (resume && !update && is.null(final_metadata) && !any(file.exists(c(outputs, processed)))) {
    candidate <- if (file.exists(final_file)) final_file else legacy_file
    if (file.exists(candidate)) {
      out <- .biohelper_adopt_legacy_final(candidate, final_file, ids, configuration, checkpoints, overwrite)
      if (!is.null(out)) return(invisible(out))
    }
  }
  if (file.exists(legacy_file)) message("Legacy final filename detected; the new result will be written to blastn_taxo_assignment.csv. Legacy file is retained.")
  executable <- .biohelper_blast_executable(blastapp_path)
  stages <- list()
  raw <- list()
  for (search in methods) {
    signature <- list(schema = "blast_search_v1", query = query_fingerprint, reference = reference,
      executable = executable, task = search, options = options[[search]], nthreads = nthreads)
    validate <- function(path) {
      x <- .biohelper_read_blast_stage(path, allow_empty = file.exists(file.path(checkpoints, paste0(search, "_search.rds"))))
      if (any(!x$qseqid %in% ids)) stop("BLAST output contains query IDs absent from the FASTA.", call. = FALSE)
      list(hit_count = nrow(x), query_count = length(unique(x$qseqid)), query_ids = sort(unique(x$qseqid)))
    }
    stage <- .biohelper_stage(paste0(search, "_search"), paste0(search, " search"), outputs[[search]], signature,
      run = function(path) .biohelper_run_assignment_search(executable$path, queries, db, search, options[[search]], nthreads,
        path, file.path(output_path, paste0(search, "_search.log"))), validate = validate,
      checkpoint_dir = checkpoints, resume = resume, overwrite = overwrite)
    raw[[search]] <- stage$value
    stages[[paste0(search, "_search")]] <- stage$metadata
  }
  # One update/validation per invocation; downstream lcaPident calls never update.
  database <- .biohelper_prepare_lca_db(update, taxonomy_db)
  tables <- list()
  for (search in methods) {
    signature <- list(schema = "lca_assignment_v2", raw = stages[[paste0(search, "_search")]]$output,
      parameters = parameters, taxonomy = database[c("identity", "provenance")], versions = versions,
      raw_provenance_verified = stages[[paste0(search, "_search")]]$provenance_verified,
      raw_generation = stages[[paste0(search, "_search")]]$generation,
      implementation = .biohelper_implementation_hash(list(lcaPident, .biohelper_lca_lookup, .biohelper_representative_order)))
    stage <- .biohelper_stage(paste0(search, "_taxonomy"), paste0(search, " taxonomy"), processed[[search]], signature,
      run = function(path) {
        result <- do.call(lcaPident, c(list(blast_file = outputs[[search]], output_file = path,
          taxonomy_db = database$path, update = FALSE, verbose = TRUE), parameters))
        provenance <- attr(result, "provenance")
        if (!identical(provenance$database$identity, database$identity)) stop("Taxonomy database changed during this workflow; rerun to resume against one database version.", call. = FALSE)
        provenance
      }, validate = function(path) {
        result <- .biohelper_read_assignment_csv(path, processed = TRUE, ids = ids)
        if (!identical(sort(result$ASVs), raw[[search]]$query_ids)) stop("Processed table is incomplete: ASV IDs differ from raw BLAST queries.", call. = FALSE)
        result
      },
      checkpoint_dir = checkpoints, resume = resume, overwrite = overwrite)
    tables[[search]] <- stage$value
    stages[[paste0(search, "_taxonomy")]] <- stage$metadata
  }
  signature <- list(schema = "blast_consensus_v2", inputs = lapply(methods, function(x) stages[[paste0(x, "_taxonomy")]]$output),
    input_generations = lapply(methods, function(x) stages[[paste0(x, "_taxonomy")]]$generation),
    queries = query_fingerprint, version = versions$biohelper,
    implementation = .biohelper_implementation_hash(list(.biohelper_merge_blast_taxonomy_tables,
      .biohelper_standardise_taxonomy_columns, .biohelper_format_blastn_taxo_assignment_output)))
  stage <- .biohelper_stage("consensus", "Consensus merging", final_file, signature,
    run = function(path) {
      merged <- .biohelper_merge_blast_taxonomy_tables(tables)
      out <- .biohelper_pad_assignments(.biohelper_format_blastn_taxo_assignment_output(merged), ids)
      readr::write_csv(out, path, na = "NA")
      list(features = nrow(out), resolution = table(rowSums(!is.na(out[, tolower(.biohelper_taxonomy_output_ranks()), drop = FALSE]))),
           assignment_methods = table(out$assignment_method, useNA = "ifany"), support_methods = table(out$support_method, useNA = "ifany"))
    }, validate = function(path) .biohelper_read_assignment_csv(path, ids = ids), checkpoint_dir = checkpoints,
    resume = resume, overwrite = overwrite, adopt = FALSE)
  stages$consensus <- stage$metadata
  provenance <- list(complete = TRUE, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    versions = versions, blast = executable, reference_database = reference,
    taxonomy_database = database[c("identity", "provenance")], parameters = parameters,
    searches = options[methods], query_count = length(ids), hit_counts = lapply(raw, function(x) x$hit_count),
    input_provenance_verified = all(vapply(stages, function(x) isTRUE(x$provenance_verified), logical(1))), stages = stages)
  .biohelper_atomic_rds(provenance, file.path(output_path, "blastn_taxo_assignment_provenance.rds"))
  invisible(stage$value)
}
