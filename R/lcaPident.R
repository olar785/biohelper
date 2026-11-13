#' Assign taxonomy from BLAST hits using LCA and identity thresholds
#'
#' @description
#' Process raw 12-column BLAST hits against a validated local NCBI taxonomy
#' database. Retains representative alignment statistics separately from the
#' multi-hit taxonomic consensus. No online taxonomy lookups are performed.
#'
#' @param blast_file Headerless 12-column BLAST TSV or CSV, in this order:
#' qseqid, qlen, pident, sseqid, sgi, evalue, bitscore, staxids, sscinames,
#' sblastnames, qcovs, qcovhsp. Empty successful BLAST output is accepted.
#' @param output_file Genuine CSV output path. Existing output is preserved
#' in an archive before a validated replacement is published.
#' @param ftbl_file Retained for compatibility; currently unused. No feature
#' table join is performed. The workflow wrapper retains every FASTA query.
#' @param minSim Species identity threshold, default 97 percent. Does not
#' filter all hits by identity. It is used only with pident=before/after.
#' @param minCov Minimum qcovs (query coverage), default 80 percent.
#' @param pident `"no"` (standalone default), `"before"` or `"after"`.
#' before masks each hit's rank according to that hit's identity, then computes
#' LCA. after computes LCA, then masks ranks using the representative hit's
#' identity. no applies no identity thresholds. The wrapper's default is before.
#' @param pgenus Genus identity threshold, default 95 percent.
#' @param pfamily Family identity threshold, default 87 percent.
#' @param porder Order identity threshold, default 83 percent.
#' @param pclass Class identity threshold, default 81 percent.
#' @param pphylum Phylum identity threshold, default 79 percent.
#' @param pkingdom Kingdom identity threshold, default 71 percent.
#' @param taxonly Retained for compatibility. Both settings return the processed
#' taxonomy with representative alignment statistics; statistics are not removed.
#' @param update Logical, default FALSE. TRUE performs one safe HTTPS taxonomy
#' update before this standalone call. FALSE never downloads. See
#' [blastn_taxo_assignment()] for cache protection and update details.
#' @param verbose Logical; print resolved/unresolved TaxID counts and output summary.
#' @param taxonomy_db Optional NCBI taxonomy SQLite path; accepts taxizedb or
#' taxonomizr names/nodes schemas. NULL uses an existing biohelper cache, otherwise
#' taxizedb's existing NCBI.sql. Missing/inaccessible/corrupt databases fail early.
#'
#' @details
#' Hits labelled uncultured, unidentified or environmental sample are removed,
#' and qcovs must meet minCov. All TaxIDs in a merged reference hit (semicolon or
#' comma separated) contribute lineage evidence; the original implementation
#' used only the first. Unresolved IDs contribute missing taxonomy, not a spurious
#' conflicting name. Unique IDs and shared ancestors are queried in batches with
#' a single read-only SQLite connection, rather than repeated classification calls.
#' An unresolved or identity-masked representative does not suppress known consensus from other hits.
#'
#' Consensus checks ranks from broad to narrow, retaining the unique known value
#' until the first disagreement and clearing that rank and all finer ranks.
#' This corrects the former fine-to-broad loop, which could retain a conflicting
#' genus or kingdom after encountering a species disagreement. The existing
#' per-rank threshold methodology, treatment of missing rank values and numeric
#' defaults are otherwise preserved.
#'
#' The representative is selected among qualifying hits by lowest E-value,
#' highest bitscore, highest identity, highest qcovs and subject identifier, with
#' a canonical row tie-breaker for remaining ties. All its statistics come from
#' the same hit, independently of the multi-hit consensus. Deterministic tie-breaking
#' can also affect after-mode rank masking when tied E-values have different
#' identities. These alignment statistics are not LCA confidence measurements.
#'
#' ASVs present in the raw table but without qualifying hits are retained with
#' missing assignment/statistics. A wholly empty file returns a typed zero-row
#' table. ASVs absent from raw output cannot be reconstructed by this standalone
#' function; the workflow wrapper adds them from its query FASTA.
#'
#' @return A data.frame written as genuine CSV. Legacy columns `ASVs`, `taxonomy`,
#' `Percent_Identity`, `Sequence_coverage` and `staxid` are retained. Added numeric
#' columns are `pident`, `qcovs`, `evalue` and `bitscore`, plus representative
#' `sseqid`. Percent_Identity/Sequence_coverage alias pident/qcovs. staxid retains
#' every TaxID attached to the representative hit. taxonomy keeps the compatible
#' nine-slot semicolon format (domain, superkingdom, kingdom, phylum, class,
#' order, family, genus, species); modern domain and legacy superkingdom are
#' normalised to domain, with missing slots written as NA.
#'
#' Attribute `provenance` records database identity/provenance, thresholds,
#' raw/qualifying hit counts, feature count and resolved/unresolved unique TaxID
#' counts. The wrapper records the method-specific assignment/support provenance
#' in its final rank table. This function always returns the table, including
#' when taxonly=FALSE, preserving the historical return type.
#'
#' @export
#' @examples
#' \dontrun{
#' result <- lcaPident(
#'   blast_file = "blastn_output.csv", output_file = "blastn_output_processed.csv",
#'   taxonomy_db = "accessionTaxa.sql", pident = "before"
#' )
#' attr(result, "provenance")
#' }
#'

lcaPident <- function(
    blast_file, output_file, ftbl_file = NULL,
    minSim = 97, minCov = 80, pident = "no",
    pgenus = 95, pfamily = 87, porder = 83,
    pclass = 81, pphylum = 79, pkingdom = 71,
    taxonly = TRUE, update = FALSE, verbose = FALSE,
    taxonomy_db = NULL
) {
  pident <- match.arg(pident, c("before", "after", "no"))
  for (path in list(blast_file, output_file)) if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) stop("blast_file and output_file must be single file paths.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.", call. = FALSE)
  thresholds <- c(kingdom = pkingdom, phylum = pphylum, class = pclass,
                  order = porder, family = pfamily, genus = pgenus, species = minSim)
  for (value in c(as.list(thresholds), list(minCov))) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value > 100) {
      stop("Identity and coverage thresholds must be numbers between 0 and 100.", call. = FALSE)
    }
  }
  database <- .biohelper_prepare_lca_db(update, taxonomy_db)
  blast <- .biohelper_read_blast_stage(blast_file)
  query_ids <- sort(unique(blast$qseqid))
  ranks <- c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  keep <- !grepl("uncultured|unidentified|environmental sample", blast$sscinames, ignore.case = TRUE) & blast$qcovs >= minCov
  keep[is.na(keep)] <- FALSE
  filtered <- blast[keep, , drop = FALSE]
  # A merged BLAST reference can represent several organisms. Use every TaxID
  # for consensus, while retaining one set of alignment statistics per hit.
  taxids <- lapply(filtered$staxids, function(value) {
    if (is.na(value)) return(character())
    unique(trimws(strsplit(as.character(value), "[;,]", perl = TRUE)[[1L]]))
  })
  ids <- unique(unlist(taxids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  lookup <- .biohelper_lca_lookup(ids, database, verbose)
  result <- data.frame(ASVs = query_ids, taxonomy = rep("Unknown", length(query_ids)),
    pident = rep(NA_real_, length(query_ids)), qcovs = rep(NA_real_, length(query_ids)),
    evalue = rep(NA_real_, length(query_ids)), bitscore = rep(NA_real_, length(query_ids)),
    sseqid = rep(NA_character_, length(query_ids)), staxid = rep(NA_character_, length(query_ids)),
    stringsAsFactors = FALSE)
  grouped <- split(seq_len(nrow(filtered)), filtered$qseqid)
  for (query in names(grouped)) {
    indexes <- grouped[[query]]
    group <- filtered[indexes, , drop = FALSE]
    best <- group[.biohelper_representative_order(group)[[1L]], , drop = FALSE]
    row <- match(query, query_ids)
    for (column in c("pident", "qcovs", "evalue", "bitscore", "sseqid")) result[[column]][[row]] <- best[[column]][[1L]]
    result$staxid[[row]] <- best$staxids[[1L]]
    evidence <- lapply(seq_along(indexes), function(i) {
      id <- taxids[[indexes[[i]]]]
      if (!length(id)) return(matrix(NA_character_, 1L, length(ranks), dimnames = list(NULL, ranks)))
      lineages <- lookup$taxonomy[id, ranks, drop = FALSE]
      if (pident == "before") {
        for (rank in names(thresholds)) if (group$pident[[i]] < thresholds[[rank]]) lineages[, rank] <- NA_character_
      }
      lineages
    })
    evidence <- do.call(rbind, evidence)
    assignment <- stats::setNames(rep(NA_character_, length(ranks)), ranks)
    # Check broad to narrow: a finer disagreement must never conceal a broader
    # conflict. Missing/unresolved lineages do not invent a conflicting name.
    for (rank in ranks) {
      values <- unique(evidence[, rank])
      values <- values[!is.na(values) & nzchar(values) & values != "NA"]
      if (length(values) > 1L) break
      if (length(values) == 1L) assignment[[rank]] <- values[[1L]]
    }
    if (pident == "after") {
      for (rank in names(thresholds)) if (best$pident[[1L]] < thresholds[[rank]]) assignment[[rank]] <- NA_character_
    }
    result$taxonomy[[row]] <- paste(ifelse(is.na(assignment), "NA", assignment), collapse = ";")
  }
  # Retain the legacy aliases for direct lcaPident callers, alongside consistent
  # numeric names used by the workflow and final output.
  result$Percent_Identity <- result$pident
  result$Sequence_coverage <- result$qcovs
  result <- result[, c("ASVs", "taxonomy", "Percent_Identity", "Sequence_coverage", "staxid",
                       "pident", "qcovs", "evalue", "bitscore", "sseqid")]
  attr(result, "provenance") <- list(database = database[c("identity", "provenance")],
    parameters = list(minSim = minSim, minCov = minCov, pident = pident, thresholds = thresholds),
    hit_count = nrow(blast), qualifying_hit_count = nrow(filtered), feature_count = length(query_ids),
    resolved_taxids = lookup$resolved, unresolved_taxids = lookup$unresolved)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".lca_", tmpdir = dirname(output_file))
  on.exit(unlink(temporary))
  readr::write_csv(result, temporary, na = "NA")
  .biohelper_atomic_publish(temporary, output_file, overwrite = FALSE)
  if (verbose) message("Taxonomy written: ", output_file, " (", length(query_ids), " ASVs)")
  result
}
