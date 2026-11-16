#' Normalise ordered taxonomy against NCBI
#'
#' @description
#' Normalise existing taxonomic names using a local NCBI taxonomy SQLite
#' database. Input names must run from broad to narrow. Their column labels
#' and number do not need to correspond to NCBI ranks. Optionally use ASV
#' sequences and BLAST to distinguish competing NCBI meanings of a name.
#'
#' @param obj Taxonomy matrix, data.frame or phyloseq object. Supply names in
#' broad-to-narrow column order, with any number of columns. Feature IDs can be
#' row names or an identifier column such as `feature_id`, `ASV` or `OTU`.
#' Semicolon-separated lineages in a `taxon` or `taxonomy` column are also accepted.
#' Recognised BLAST alignment/provenance columns are preserved separately and
#' excluded from the ordered-name search; they are not taxonomy ranks.
#' @param sqlFile Path to the local NCBI taxonomy SQLite database, containing
#' names and nodes (for example `accessionTaxa.sql`). This is used for name lookup,
#' ancestry and output ranks. It is separate from the BLAST sequence database `db`.
#' @param addExtra Logical, default TRUE. Apply the existing extra kingdom-group
#' assignments, including Protozoa and Archaeplastida (excluding Viridiplantae),
#' where applicable. Set FALSE to retain the NCBI kingdom assignments directly.
#' @param spnc Logical, default FALSE. If TRUE, join the final two taxonomic
#' columns as genus and species epithet. This is an explicit formatting request,
#' independent of column labels; it requires at least two taxonomic columns.
#' Otherwise supply complete binomial names where available.
#' @param ranks NCBI ranks to return. Default: Domain, Kingdom, Phylum, Class,
#' Order, Family, Genus and Species. Output rank column names are lowercase.
#' This controls output only; it does not constrain the input names or matching.
#' @param sequences ASV sequences as a FASTA file path, named DNAStringSet or
#' uniquely named character vector. IDs must match the taxonomy feature IDs;
#' FASTA headers are read up to the first whitespace. For phyloseq input,
#' reference sequences are used automatically when available. ASVs without
#' matching sequences are recorded as `missing_sequence` and use the configured
#' ambiguity fallback. You do not need to intersect IDs manually.
#' @param blast_hits Optional precomputed raw BLAST data.frame, or 12-column
#' file readable by [read_blast_taxo_csv()]. Required columns: `qseqid`, `sseqid`,
#' `staxids`, `pident`, `qcovhsp`, `evalue` and `bitscore`. No search is run when
#' supplied, so `sequences`, `db` and `blastapp_path` are not required.
#' @param blast_options Named list of optional BLAST settings. Omitted entries
#' use these defaults:
#' \describe{
#'   \item{`scope = "ambiguous"`}{Which assignments to check: `"ambiguous"`
#'     checks names still matching multiple NCBI candidates; `"all"` also audits
#'     names with one candidate. Neither option searches features without a
#'     recognised name candidate.}
#'   \item{`fallback = "first"`}{For remaining ambiguous names: `"first"`,
#'     `"lca"` or `"unassigned"`. Use `"lca"` for BLAST followed by LCA.}
#'   \item{`remote = FALSE`}{Run BLAST locally by default. TRUE submits selected
#'     DNA sequences to NCBI's online BLAST service; if `db` is omitted, that
#'     service uses `core_nt`. No local executable is needed for online searches.}
#'   \item{`task = "blastn"`}{BLAST task: `"blastn"` or `"megablast"`.}
#'   \item{`nthreads = 1`}{Number of threads for local searches.}
#'   \item{`verbose = TRUE`}{Show name-matching counts, sequence deduplication,
#'     cache reuse, search settings and completion/fallback summaries. FALSE
#'     suppresses progress messages, but does not suppress warnings.}
#'   \item{`progress_seconds = 30`}{Interval in seconds between elapsed-time
#'     updates while local BLAST is running. BLAST does not report a reliable
#'     percentage or completion-time estimate.}
#'   \item{`max_hits = 100`}{Maximum reference hits retained per query.}
#'   \item{`evalue = 1e-10`}{Maximum E-value for search and hit screening.}
#'   \item{`min_identity = 95`}{Minimum alignment identity, in percent.}
#'   \item{`min_coverage = 90`}{Minimum query coverage, in percent, of one
#'     alignment segment (HSP), rather than pooled short alignments.}
#'   \item{`score_fraction = 0.98`}{Retain hits scoring at least this fraction
#'     of the best hit passing identity, coverage and E-value screening.}
#'   \item{`min_accessions = 2`}{Minimum number of distinct reference accessions.}
#'   \item{`min_taxa = 2`}{Minimum number of distinct reference TaxIDs among
#'     retained hits. Setting 1 allows evidence from a single TaxID while still
#'     requiring `min_accessions`.}
#'   \item{`min_support = 1`}{Required fraction of retained hit TaxIDs supporting
#'     a candidate lineage. Default 1 requires all; must be greater than 0.5.}
#'   \item{`conflict_action = "flag"`}{When qualifying hits conflict with all
#'     name candidates, flag for review and apply the ambiguity fallback.
#'     `"unassign"` explicitly clears the assignment instead. Unique audited
#'     names remain assigned under `"flag"`; they have no competing-name LCA.}
#'   \item{`feature_ids = NULL`}{Optional feature subset. By default all eligible
#'     features are considered and sequence availability is matched automatically.}
#'   \item{`cache_dir = file.path(tempdir(), "biohelper-blast")`}{Search cache.
#'     Use a persistent directory to reuse raw hits or resume online requests
#'     between R sessions. Identical DNA reuses hits; different sequences are
#'     never pooled merely because their taxonomy names agree.}
#'   \item{`refresh = FALSE`}{Reuse matching cached searches; TRUE reruns them.}
#'   \item{`max_remote_queries = 25`}{Maximum uncached distinct DNA sequences
#'     per online call. Does not limit local searches.}
#'   \item{`poll_seconds = 60`}{Minimum interval between online job checks.}
#'   \item{`timeout_seconds = 900`}{Maximum online waiting time per call.}
#'   \item{`email = NULL`}{Optional contact email for online requests.}
#' }
#' @param return_details Logical, default FALSE. FALSE returns the usual taxonomy
#' table or phyloseq object. TRUE changes the return value to a list containing
#' `taxonomy`, `blast` and `ambiguity`; see Value for the complete output contract.
#' @param ambiguity How to handle competing names after ordered lineage matching:
#' `"first"` (default), `"lca"`, `"unassigned"` or `"blast_check"`.
#' Only `"blast_check"` enables sequence checking. The other strategies do not
#' run BLAST, access its sequence database or submit sequences online.
#' @param db BLAST sequence database, using the same convention as
#' [blastn_taxo_assignment()]. For local searches, supply the formatted database
#' prefix, e.g. `/databases/nt`, not a `.nin` file, FASTA file or SQLite file.
#' A database name resolvable through BLAST's database search path also works.
#' Reference records must have TaxIDs for lineage checking. For online searches,
#' supply an NCBI database name, or omit it to use `core_nt`. Default NULL;
#' required only when a new local search actually runs.
#' @param blastapp_path Path to the `blastn` executable, using the same name as
#' [blastn_taxo_assignment()]. Default: detect `blastn` on the system PATH.
#' Required only for a new local search; ignored for online/precomputed evidence.
#'
#' @details
#' Names are searched right-to-left, beginning with the finest supplied name.
#' Unrecognised names fall back to broader names. For a name with competing
#' NCBI IDs, broader supplied names are compared against the complete NCBI
#' lineages, including clades and intermediate ranks. Supporting matches must
#' follow the input ancestry order; repeated names can match the same node.
#' Names unable to support an ordered match in any remaining candidate are skipped.
#' Matching stops when one candidate remains. NCBI supplies the output ranks.
#' Identical ordered taxonomy profiles share name/lineage matching work, but
#' each feature retains its own sequence evidence and diagnostics. All selected
#' sequences are collected before BLAST; identical DNA is searched once and
#' previously cached searches are reused. New local searches run in one batch.
#'
#' With `ambiguity = "first"`, remaining supported candidates use the first
#' lookup ID with a warning. If no broader name supports any candidate, this
#' default policy leaves the assignment unassigned. `"lca"` instead uses the
#' deepest common NCBI ancestor of the remaining candidates; `"unassigned"`
#' clears the ambiguous assignment. Requested ranks below an LCA remain NA.
#' LCA nodes can be clades not present among the requested output ranks.
#'
#' With `ambiguity = "blast_check"`, sufficiently similar, long, near-best
#' sequence hits are checked against the candidate TaxIDs. Each distinct hit
#' TaxID contributes one vote. The deepest supported candidate is selected when
#' candidates are nested. This selects a meaning of the supplied name; it does
#' not automatically replace it with the species name of a BLAST hit.
#' Weak, mixed, missing or unavailable evidence uses `blast_options$fallback`.
#' The screening thresholds are adjustable, not validated rank-specific cutoffs;
#' multiple reference accessions need not be independent biological observations.
#'
#' A failed local search, missing database/executable, network error or pending
#' online job stops the call; these execution errors are not treated as LCA
#' evidence. Pending online request IDs are saved. Rerun the same call with the
#' same persistent cache to resume. Completed hits are cached separately from
#' screening thresholds, allowing thresholds to be changed without a new search.
#'
#' @return
#' With `return_details = FALSE` (default): matrix/data.frame input returns a
#' data.frame with `feature_id` followed by requested lowercase rank columns;
#' phyloseq input returns a phyloseq object with normalised taxonomy. The input
#' object is not modified. Unassigned ranks are NA; an unassigned requested
#' domain is represented as `"Unknown"`.
#'
#' When BLAST checking is enabled, attribute `blast_check` contains its report.
#' With a non-default ambiguity strategy, attribute `ambiguity` contains the
#' remaining-candidate/fallback diagnostics. These are attributes, not extra
#' taxonomy columns. Ordinary `"first"` calls retain the original output shape.
#'
#' With `return_details = TRUE`, the return value is explicitly a list:
#' \describe{
#'   \item{`taxonomy`}{The normalised data.frame or phyloseq object.}
#'   \item{`blast`}{NULL when BLAST is disabled; otherwise a list containing
#'     `diagnostics` (one row per selected feature), raw `hits`, applied
#'     `settings`, and `hit_source` (`"blastn"` or `"precomputed"`). Diagnostics
#'     include feature/candidate/selected TaxIDs, final ID, status, hit counts,
#'     support fraction, minimum retained identity and coverage. Status examples:
#'     `blast_resolved`, `blast_supported`, `missing_sequence`, `no_hits`,
#'     `insufficient_quality`, `insufficient_records`, `mixed_hits`,
#'     `taxonomy_unavailable`, `insufficient_resolution` and `blast_conflict`.}
#'   \item{`ambiguity`}{One row per unresolved competing-name feature handled by
#'     a fallback, with `feature_id`, `candidate_ids`, `action`, `selected_id`,
#'     `lca_id`, `lca_name`, `lca_rank` and requested `method`. An empty table
#'     means no competing names remain for that fallback.}
#' }
#' Ambiguity warnings report a feature count and up to five example IDs; their
#' `feature_ids` field retains all affected IDs. BLAST review warnings identify
#' features whose evidence did not resolve or confirm the name assignment.
#'
#' @export
#' @examples
#' \dontrun{
#' # Name-based normalisation: no BLAST required.
#' normalised <- taxo_normalisation(taxa, sqlFile = "accessionTaxa.sql")
#'
#' # Local BLAST with LCA fallback. Other BLAST options use their defaults.
#' normalised <- taxo_normalisation(
#'   taxa, sqlFile = "accessionTaxa.sql", sequences = "uniqueSeqs.fasta",
#'   ambiguity = "blast_check", db = "/databases/nt",
#'   blast_options = list(fallback = "lca")
#' )
#'
#' # Explicitly request the additional reports.
#' details <- taxo_normalisation(
#'   taxa, sqlFile = "accessionTaxa.sql", sequences = "uniqueSeqs.fasta",
#'   ambiguity = "blast_check", db = "/databases/nt",
#'   blast_options = list(fallback = "lca"), return_details = TRUE
#' )
#' details$taxonomy
#' details$blast$diagnostics
#' details$ambiguity
#'
#' # Online BLAST is opt-in; db defaults to core_nt for this mode only.
#' details <- taxo_normalisation(
#'   taxa, sqlFile = "accessionTaxa.sql", sequences = "uniqueSeqs.fasta",
#'   ambiguity = "blast_check",
#'   blast_options = list(remote = TRUE, fallback = "lca"), return_details = TRUE
#' )
#' }
#'

taxo_normalisation = function(obj, sqlFile, addExtra = TRUE, spnc = FALSE, ranks = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                              sequences = NULL,
                              blast_hits = NULL, blast_options = list(), return_details = FALSE,
                              ambiguity = c("first", "lca", "unassigned", "blast_check"),
                              db = NULL, blastapp_path = unname(Sys.which("blastn"))){
  ambiguity <- match.arg(ambiguity)
  use_blast <- ambiguity == "blast_check"
  if (!is.logical(return_details) || length(return_details) != 1L || is.na(return_details)) {
    stop("`return_details` must be TRUE or FALSE.", call. = FALSE)
  }
  if (use_blast) {
    blast_options <- .biohelper_blast_options(blast_options)
    if (is.null(sequences) && is.null(blast_hits) && inherits(obj, "phyloseq")) {
      sequences <- phyloseq::refseq(obj, errorIfNULL = FALSE)
    }
  }

  prepared <- .biohelper_prepare_normalisation(
    obj, sqlFile, spnc, ranks, progress = if (use_blast) blast_options else NULL
  )
  .biohelper_finish_normalisation(
    prepared, sqlFile, addExtra, ambiguity, sequences, blast_hits, blast_options,
    return_details, db, blastapp_path
  )
}

.biohelper_prepare_normalisation <- function(obj, sqlFile, spnc, ranks,
                                             cache = .biohelper_normalisation_cache(), progress = NULL) {
  `%ni%` <- Negate(`%in%`)

  if("phyloseq" %in% class(obj)){
    df = obj@tax_table@.Data %>% as.data.frame()
  } else{
    df = obj %>% as.data.frame()
  }

  ranks = tolower(ranks)
  input_labels <- colnames(df)
  missing_labels <- is.na(input_labels) | input_labels == ""
  input_labels[missing_labels] <- paste0("input_taxon_", which(missing_labels))
  colnames(df) <- make.unique(tolower(input_labels))
  taxa_id = c("otu","otus","asv","asvs","feature_id","feature.id","nR")
  # Alignment/provenance fields can accompany a taxonomy table, but must never
  # become names in the ordered taxonomy search. Restore them without cleaning.
  metadata_columns <- intersect(c(.biohelper_blast_stats_columns(), "bitscore", "sseqid"), colnames(df))
  feature_metadata <- df[, metadata_columns, drop = FALSE]
  id_column <- .biohelper_detect_column(df, setdiff(taxa_id, "nR"))
  metadata_ids <- if (is.na(id_column)) rownames(df) else as.character(df[[id_column]])
  df <- df[, !colnames(df) %in% metadata_columns, drop = FALSE]

  # Recognise a packed lineage by its separator, not just its column label:
  # a plain name in a wide column named taxon/taxonomy remains a normal clue.
  packed_columns <- intersect(c("taxon", "taxonomy"), colnames(df))
  packed_columns <- packed_columns[vapply(packed_columns, function(column) {
    any(grepl(";", df[[column]], fixed = TRUE), na.rm = TRUE)
  }, logical(1))]
  if (length(packed_columns)) {
    column <- packed_columns[[1L]]
    df[[column]] <- gsub("D_\\d+__", "", df[[column]])
    df <- df[, colnames(df) %in% c(taxa_id, column), drop = FALSE]
    df <- splitstackshape::cSplit(df, column, sep = ";")
  }
  df = df %>% dplyr::mutate_if(is.factor, as.character)
  df[is.na(df)]<-""

  if(any(taxa_id %in% colnames(df))){
    taxa_id_df = which(colnames(df) %in% taxa_id)%>% as.numeric()
    colnames(df)[taxa_id_df] = "feature_id"
    if(any(duplicated(df))){
      cat("Duplicated rows! Keeping the first instance only.\n")
      df = df[!duplicated(df), ]
      rownames( df ) <- NULL
      df = df %>% column_to_rownames("feature_id")
    }else{
      rownames( df ) <- NULL
      df = df %>% column_to_rownames("feature_id")
    }
  }

  paternsToRemove = c("^.+_environmental.+|environmental_.+|uncultured_.+|_sp\\..+|_sp.|_sp.+| sp\\..+| sp.| sp.+|_\\(.+|^.+_metagenome|_cf.|s__|g__|f__|c__|o__|c__|p__|'")
  df = df %>% dplyr::mutate_all(list(~str_replace(.,paternsToRemove, ""))) %>% dplyr::mutate_all(list(~na_if(.,"")))
  df = df %>%
    dplyr::mutate_all(list(~str_replace(.,"NA NA| NA", ""))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), function(x) gsub(pattern = "_", replacement = " ", x = x))) %>%
    dplyr::mutate_all(list(~na_if(.,""))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), function(x) gsub(pattern = "^\\s+|\\s+$", replacement = "", x = x))) # The latter is to removing leading and trailing spaces

  ranks_indexes = which(colnames(df) %ni% c("otu","otus","asv","asvs","feature_id","feature.id","nR"))
  non_taxo_ranks = c("otu","otus","asv","asvs","feature_id","feature.id","nR")
  if (nrow(df) == 0L || length(ranks_indexes) == 0L) {
    stop("`obj` must contain at least one feature and one taxonomic column.", call. = FALSE)
  }
  if (!is.logical(spnc) || length(spnc) != 1L || is.na(spnc)) {
    stop("`spnc` must be TRUE or FALSE.", call. = FALSE)
  }
  if (spnc) {
    if (length(ranks_indexes) < 2L) {
      stop("`spnc = TRUE` requires at least two taxonomic columns.", call. = FALSE)
    }
    last <- utils::tail(ranks_indexes, 2L)
    df[[last[[2L]]]] <- trimws(paste(
      ifelse(is.na(df[[last[[1L]]]]), "", df[[last[[1L]]]]),
      ifelse(is.na(df[[last[[2L]]]]), "", df[[last[[2L]]]])
    ))
  }
  original_df <- df
  profile_keys <- .biohelper_taxonomy_profile_keys(df[, ranks_indexes, drop = FALSE])
  representatives <- which(!duplicated(profile_keys))
  profile_map <- match(profile_keys, profile_keys[representatives])
  df <- df[representatives, , drop = FALSE]
  .biohelper_blast_message(progress, "NCBI name matching: ", nrow(original_df),
    " features; ", nrow(df), " distinct ordered taxonomy profiles.")
  # Input columns express the intended resolution order, not authoritative NCBI
  # ranks. Search non-missing names from right to left and retain the position
  # of the name that actually matches, including after an unsuccessful lookup.
  taxonomy <- as.matrix(df[, ranks_indexes, drop = FALSE])
  available <- !is.na(taxonomy) & nzchar(taxonomy)
  lookup_positions <- lapply(seq_len(nrow(df)), function(row) {
    rev(ranks_indexes[available[row, ]])
  })
  lookup_step <- rep(1L, nrow(df))
  res_df <- data.frame(
    feature_id = rownames(df), rpt_indexes = NA_integer_,
    taxa = NA_character_, id = NA_character_, stringsAsFactors = FALSE
  )
  pending <- which(lengths(lookup_positions) > 0L)
  while (length(pending) > 0L) {
    positions <- vapply(pending, function(row) {
      lookup_positions[[row]][[lookup_step[[row]]]]
    }, integer(1))
    taxa <- vapply(seq_along(pending), function(j) {
      df[[positions[[j]]]][[pending[[j]]]]
    }, character(1))
    res_df$rpt_indexes[pending] <- positions
    res_df$taxa[pending] <- taxa
    .biohelper_blast_message(progress, "Looking up ", length(unique(taxa)),
      " distinct names for ", length(pending), " unresolved profiles.")
    res_df$id[pending] <- .biohelper_cached_ncbi_names(taxa, sqlFile, cache)
    pending <- pending[is.na(res_df$id[pending])]
    lookup_step[pending] <- lookup_step[pending] + 1L
    pending <- pending[lookup_step[pending] <= lengths(lookup_positions[pending])]
  }
  # Resolve ambiguous names using broader names in positional order, without
  # interpreting column labels as ranks or limiting the number of columns.
  ambiguous_rows <- which(!is.na(res_df$id) & grepl(",", res_df$id, fixed = TRUE))
  blast_candidates <- lapply(res_df$id, function(id) {
    if (is.na(id)) character() else unique(trimws(strsplit(id, ",", fixed = TRUE)[[1]]))
  })
  ambiguous_features <- character()
  if (length(ambiguous_rows) > 0L) {
    candidate_ids <- lapply(
      strsplit(res_df$id[ambiguous_rows], ",", fixed = TRUE),
      function(ids) unique(trimws(ids))
    )
    lineage_ids <- unique(unlist(candidate_ids))
    # getRawTaxonomy formats numeric IDs and can pad shorter names with spaces.
    # Its results follow the requested order, so name them with our exact IDs.
    .biohelper_blast_message(progress, "Checking full NCBI lineages for ",
      length(ambiguous_rows), " ambiguous profiles (", length(lineage_ids), " candidate TaxIDs).")
    candidate_lineages <- .biohelper_cached_ncbi_lineages(lineage_ids, sqlFile, cache)
    ambiguous_features <- character()
    for (j in seq_along(ambiguous_rows)) {
      row <- ambiguous_rows[[j]]
      candidates <- candidate_ids[[j]]
      has_match <- FALSE
      # Raw lineages run from the candidate towards the root. Each successive
      # input clue must match at or above that candidate's previous clue.
      lineage_positions <- stats::setNames(rep(1L, length(candidates)), candidates)
      higher_columns <- rev(ranks_indexes[ranks_indexes < res_df$rpt_indexes[[row]]])
      for (column in higher_columns) {
        supplied_taxon <- df[[column]][[row]]
        if (is.na(supplied_taxon) || supplied_taxon == "") {
          next
        }
        positions <- vapply(candidates, function(id) {
          found <- which(unname(candidate_lineages[[id]]) == supplied_taxon)
          found <- found[found >= lineage_positions[[id]]]
          if (length(found)) found[[1L]] else NA_integer_
        }, integer(1))
        matching <- candidates[!is.na(positions)]
        if (length(matching) > 0L) {
          has_match <- TRUE
          lineage_positions[matching] <- positions[!is.na(positions)]
          candidates <- matching
          if (length(candidates) == 1L) {
            break
          }
        }
      }
      blast_candidates[[row]] <- candidates
      res_df$id[[row]] <- if (has_match) candidates[[1L]] else NA_character_
      if (has_match && length(candidates) > 1L) {
        ambiguous_features <- c(ambiguous_features, res_df$feature_id[[row]])
      }
    }
  }
  ambiguous_profiles <- res_df$feature_id %in% ambiguous_features
  res_df <- res_df[profile_map, , drop = FALSE]
  rownames(res_df) <- NULL
  res_df$feature_id <- rownames(original_df)
  blast_candidates <- blast_candidates[profile_map]
  ambiguous_features <- res_df$feature_id[ambiguous_profiles[profile_map]]
  .biohelper_blast_message(progress, "Name matching complete: ",
    sum(lengths(blast_candidates) == 1L), " unique matches; ",
    sum(lengths(blast_candidates) > 1L), " ambiguous ASVs; ",
    sum(lengths(blast_candidates) == 0L), " without a recognised name.")
  list(obj = obj, df = original_df, res_df = res_df, blast_candidates = blast_candidates,
    ambiguous_features = ambiguous_features, ranks = ranks, non_taxo_ranks = non_taxo_ranks,
    metadata_columns = metadata_columns, feature_metadata = feature_metadata, metadata_ids = metadata_ids)
}

.biohelper_finish_normalisation <- function(prepared, sqlFile, addExtra, ambiguity,
                                             sequences, blast_hits, blast_options,
                                             return_details, db, blastapp_path,
                                             blast_context = NULL) {
  obj <- prepared$obj
  df <- prepared$df
  res_df <- prepared$res_df
  blast_candidates <- prepared$blast_candidates
  ambiguous_features <- prepared$ambiguous_features
  ranks <- prepared$ranks
  non_taxo_ranks <- prepared$non_taxo_ranks
  metadata_columns <- prepared$metadata_columns
  feature_metadata <- prepared$feature_metadata
  metadata_ids <- prepared$metadata_ids
  use_blast <- ambiguity == "blast_check"
  fallback_action <- if (use_blast) blast_options$fallback else ambiguity
  blast_result <- NULL
  resolved <- character()
  blocked_by_blast <- character()
  if (use_blast) {
    blast_result <- .biohelper_check_blast(
      res_df, blast_candidates, sqlFile, blast_options$scope, sequences,
      if (is.null(blast_context)) blast_hits else blast_context$hits, blast_options,
      db = db, blastapp_path = blastapp_path,
      missing_sequences = blast_context$missing, paths = blast_context$paths
    )
    if (!is.null(blast_context)) blast_result$hit_source <- "blastn"
    res_df$id <- blast_result$ids
    diagnostics <- blast_result$diagnostics
    resolved <- diagnostics$feature_id[diagnostics$status == "blast_resolved"]
    if (blast_options$conflict_action == "unassign") {
      blocked_by_blast <- diagnostics$feature_id[diagnostics$status == "blast_conflict"]
    }
    ambiguous_features <- setdiff(ambiguous_features, resolved)
    ambiguous_features <- ambiguous_features[!is.na(res_df$id[match(ambiguous_features, res_df$feature_id)])]
    needs_review <- diagnostics$feature_id[!diagnostics$status %in% c("blast_supported", "blast_resolved")]
    if (length(needs_review)) {
      warning(structure(list(
        message = paste0("BLAST did not resolve or confirm ", length(needs_review),
          " selected feature(s). Inspect the returned BLAST diagnostics and fallback; clear conflicts configured to unassign remain unassigned."),
        call = NULL, feature_ids = needs_review
      ), class = c("biohelper_blast_review", "warning", "condition")))
    }
    blast_result$ids <- NULL
  }
  ambiguity_rows <- which(lengths(blast_candidates) > 1L &
    !res_df$feature_id %in% c(resolved, blocked_by_blast))
  ambiguity_result <- .biohelper_resolve_ambiguity(
    res_df, blast_candidates, ambiguity_rows, fallback_action, sqlFile, paths = blast_context$paths
  )
  res_df$id <- ambiguity_result$ids
  ambiguity_diagnostics <- ambiguity_result$diagnostics
  ambiguity_diagnostics$method <- rep(ambiguity, nrow(ambiguity_diagnostics))
  if (!is.null(blast_result) && nrow(blast_result$diagnostics)) {
    blast_result$diagnostics$final_id <- res_df$id[match(blast_result$diagnostics$feature_id, res_df$feature_id)]
  }
  if (fallback_action == "first" && length(ambiguous_features) > 0L) {
      examples <- utils::head(ambiguous_features, 5L)
      warning(structure(
        list(
          message = paste0(
            "Multiple NCBI IDs still match the supplied broader names for ",
            length(ambiguous_features), " feature(s): ", paste(examples, collapse = ", "),
            if (length(ambiguous_features) > 5L) " (first 5 shown)" else "",
            ". Using the first matching ID; these assignments remain ambiguous."
          ),
          call = NULL,
          feature_ids = ambiguous_features
        ),
        class = c("biohelper_taxonomy_ambiguity", "warning", "condition")
      ))
  }
  if (fallback_action != "first" && length(ambiguity_rows)) {
    warning(structure(list(
      message = paste0("Competing NCBI names remain for ", length(ambiguity_rows),
        " feature(s). ", if (fallback_action == "lca") "Retained their lowest common ancestor (LCA)." else "Left their taxonomy unassigned.",
        " Inspect the ambiguity diagnostics for candidate and selected TaxIDs."),
      call = NULL, feature_ids = res_df$feature_id[ambiguity_rows], action = fallback_action
    ), class = c("biohelper_taxonomy_ambiguity", "warning", "condition")))
  }
  if (use_blast) {
    .biohelper_blast_message(blast_options, "Resolving taxonomy: ", length(resolved),
      " features resolved by BLAST; ", length(ambiguity_rows),
      " remaining name ambiguities handled with '", fallback_action, "'.")
  }
  res_df[ranks] = taxonomizr::getTaxonomy(res_df$id, sqlFile = sqlFile, desiredTaxa = str_to_lower(ranks))
  res_df = res_df[,colnames(res_df) %in% c(ranks,non_taxo_ranks)]
  if ("domain" %in% colnames(res_df)) {
    res_df$domain = res_df$domain %>% replace_na("Unknown")
  }
  res_df = res_df %>% column_to_rownames("feature_id")
  res_df=cbind(res_df,df[,colnames(df) %in% non_taxo_ranks,drop = FALSE])

  if(addExtra & ("kingdom" %in% str_to_lower(ranks))){
    res_df = biohelper::extra_taxo_assignment(obj = res_df, sqlFile = sqlFile)
  }
  if (length(metadata_columns)) {
    res_df <- cbind(res_df, feature_metadata[match(rownames(res_df), metadata_ids), , drop = FALSE])
  }

  if("phyloseq" %in% class(obj)){
    obj@tax_table@.Data = res_df %>% as.matrix()
    output <- obj
  }else{
    output <- res_df %>% rownames_to_column("feature_id")
  }
  if (return_details) return(list(taxonomy = output, blast = blast_result, ambiguity = ambiguity_diagnostics))
  if (!is.null(blast_result)) attr(output, "blast_check") <- blast_result
  if (ambiguity != "first") attr(output, "ambiguity") <- ambiguity_diagnostics
  return(output)
}
