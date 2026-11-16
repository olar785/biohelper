#' taxo_merge
#'
#' @description
#' This function takes in a list of taxonomic assignment dataframes, normalizes the
#' taxonomy using NCBI curated taxonomic database, and uses the highest
#' taxonomic resolution between the assignments if there is consensus across
#' all assigned ranks and dataframes. Otherwise, it assigns taxonomy the Last common
#' ancestor (LCA) if there is a consensus (>50%) across dataframes.
#' By using the 'priority_df' option, it is also possible to use the taxonomy of a specific dataframe (assignment) instead of using LCA.
#' `taxo_merge()` returns lowercase rank column names and standardises the
#' highest taxonomic rank to one `domain` column in the returned table.
#'
#' Importantly, this function requires the download of NCBI taxonomic database (~65 GB).
#' Additionally, the first column of the dataframe needs to be the feature id (e.g. ASV, OTU, MOTU, ZOTU identifier, etc.).
#'
#' This takes few minutes to install using the following command:
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package.
#'
#' @param
#' df_list        List of dataframes with taxonomic assignment
#' @param
#' sqlFile        Path to NCBI taxonomic reference database
#' @param
#' ranks          Ranks to return
#' @param
#' priority_df    Taxonomic dataframe to prioritize in case of discrepancy between assignments (default is NA).
#' @param
#' addExtra       Currently adds the TSAR and Archaeplastida (excl. Viridiplantae) groups under the Kingdom rank (default is TRUE). This is to make it easier to dissociate the Plantea and TSAR group from Metazoans
#' @param
#' spnc           Only needs to be applied if the Genus is not present under the Species column (default is FALSE).
#' @param ambiguity How to resolve competing NCBI name matches within each
#' input: `"first"` (default), `"lca"`, `"unassigned"`, or `"blast_check"`.
#' Passed to [taxo_normalisation()] before merging. This is separate from
#' resolving disagreements between input tables during the merge.
#' @param sequences Shared ASV sequences as a FASTA path, named DNAStringSet,
#' or named character vector. Sequence IDs are matched separately to each
#' input's feature IDs. Used only with `ambiguity = "blast_check"`.
#' The FASTA is read once for the complete call.
#' @param db Formatted BLAST sequence database prefix or resolvable database
#' name, as in [blastn_taxo_assignment()]. Separate from `sqlFile`; local
#' reference sequences must have TaxIDs. Used only for BLAST checking.
#' @param blastapp_path Path to the `blastn` executable. Defaults to detection
#' on PATH. Used only when a new local BLAST search is needed.
#' @param blast_options Named list of settings passed to [taxo_normalisation()].
#' Local searches are the default (`remote = FALSE`), with `scope = "ambiguous"`
#' and `fallback = "lca"`. Remaining ambiguous names use their common ancestor
#' when evidence is insufficient or a matching sequence is missing.
#' A shared persistent `cache_dir` reuses searches
#' for identical sequences across tables and calls. All inputs are matched to
#' NCBI first, then their selected ASVs are combined into one BLAST batch of
#' uncached, distinct DNA sequences. Even `refresh = TRUE` searches each distinct
#' sequence once per call. `feature_ids`, when supplied, selects from the union
#' of input ASV IDs. The online `max_remote_queries` limit applies to the pooled
#' search. `verbose = TRUE` reports counts and progress, including elapsed-time
#' updates in minutes every five minutes (`progress_seconds = 300`) while local
#' BLAST runs. `nthreads = 1`
#' remains the default; supply a higher thread count to use more local CPUs.
#' See [taxo_normalisation()]
#' for all options and screening defaults. Search failures stop the call.
#' @return A data.frame with `feature_id` and requested lowercase taxonomy
#' ranks. BLAST alignment and provenance columns are excluded from merging and
#' are not returned. Each input is normalised using the supplied options;
#' the merged result is not normalised again. For a non-default `ambiguity`
#' strategy, `attr(result, "normalisation")` contains one report per input
#' (retaining list names when supplied), each with `blast` and `ambiguity`
#' diagnostics from [taxo_normalisation()]. These reports describe the input
#' normalisations, rather than the later merge decisions.
#'
#' @export
#' @examples
#' \dontrun{
#' taxo_merge(
#'   df_list = list(df1, df2),
#'   sqlFile = path_to_NCBI_taxo_db,
#'   ranks = c(
#'     "Domain", "Kingdom", "Phylum", "Class",
#'     "Order", "Family", "Genus", "Species"
#'   ),
#'   priority_df = NA,
#'   addExtra = TRUE,
#'   spnc = FALSE
#' )
#'
#' # Local BLAST: LCA fallback and five-minute progress updates are defaults.
#' consensus_taxo <- taxo_merge(
#'   df_list = list(silva = taxa_silva, pr2 = taxa_pr2, blast = blastn),
#'   sqlFile = "/databases/accessionTaxa.sql", # NCBI names and ancestry
#'   sequences = "uniqueSeqs.fasta",          # ASV IDs must match each input
#'   ambiguity = "blast_check",
#'   db = "/databases/nt",                    # BLAST sequence database prefix
#'   blastapp_path = "/usr/local/bin/blastn",
#'   blast_options = list(
#'     remote = FALSE,
#'     scope = "ambiguous",
#'     nthreads = 4,                          # choose for your computer
#'     cache_dir = "blast_cache"              # reuse hits between R sessions
#'   )
#' )
#' reports <- attr(consensus_taxo, "normalisation")
#' reports$pr2$blast$diagnostics
#' reports$pr2$ambiguity
#' }


taxo_merge = function(
    df_list,
    sqlFile,
    ranks = c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    priority_df = NA,
    addExtra = T,
    spnc = F,
    ambiguity = c("first", "lca", "unassigned", "blast_check"),
    sequences = NULL,
    db = NULL,
    blastapp_path = unname(Sys.which("blastn")),
    blast_options = list())
{
  if (!is.list(df_list) || is.data.frame(df_list) || !length(df_list)) {
    stop("`df_list` must be a non-empty list of taxonomy tables.", call. = FALSE)
  }
  ambiguity <- match.arg(ambiguity)
  if (ambiguity == "blast_check") {
    blast_options <- .biohelper_blast_options(blast_options)
  }
  ranks <- .biohelper_standardise_rank_names(ranks)
  if (!length(ranks)) {
    stop("`ranks` must include at least one recognised taxonomy rank.", call. = FALSE)
  }
  normalisation_ranks <- .biohelper_taxonomy_normalisation_ranks(ranks)
  reports <- vector("list", length(df_list))
  names(reports) <- names(df_list)
  plans <- vector("list", length(df_list))
  cache <- .biohelper_normalisation_cache()
  use_blast <- ambiguity == "blast_check"

  for (i in seq_along(df_list)) {
    if (use_blast) {
      .biohelper_blast_message(blast_options, "Planning taxonomy for table ", i, "/", length(df_list),
        if (!is.null(names(df_list))) paste0(" (", names(df_list)[[i]], ")") else "", ".")
    } else cat("\nNormalizing taxonomy of object #",i,"\n")
    colnames(df_list[[i]]) = colnames(df_list[[i]]) %>% tolower()
    colnames(df_list[[i]])[1] = "feature_id"
    if (use_blast) {
      plans[[i]] <- .biohelper_prepare_normalisation(
        df_list[[i]], sqlFile, spnc, normalisation_ranks, cache, blast_options
      )
    } else {
      df_list[[i]] <- taxo_normalisation(
        obj = df_list[[i]], sqlFile = sqlFile, ranks = normalisation_ranks,
        addExtra = addExtra, spnc = spnc, ambiguity = ambiguity,
        sequences = sequences, db = db, blastapp_path = blastapp_path,
        blast_options = blast_options
      )
    }
  }
  if (use_blast) {
    pooled <- .biohelper_pool_normalisation_blast(plans, sequences, blast_options, db, blastapp_path, sqlFile)
    for (i in seq_along(df_list)) {
      .biohelper_blast_message(blast_options, "Applying BLAST evidence and fallback to table ", i, "/", length(df_list), ".")
      options_i <- blast_options
      if (!is.null(options_i$feature_ids)) {
        options_i$feature_ids <- intersect(options_i$feature_ids, plans[[i]]$res_df$feature_id)
      }
      df_list[[i]] <- .biohelper_finish_normalisation(
        plans[[i]], sqlFile, addExtra, ambiguity, NULL, NULL, options_i,
        FALSE, db, blastapp_path, blast_context = pooled
      )
    }
  }
  for (i in seq_along(df_list)) {
    reports[[i]] <- list(
      blast = attr(df_list[[i]], "blast_check"),
      ambiguity = attr(df_list[[i]], "ambiguity")
    )
    df_list[[i]] = .biohelper_standardise_taxonomy_columns(
      df_list[[i]],
      id_col = "feature_id",
      output_id_col = "feature_id"
    ) %>%
      dplyr::mutate(df=as.character(i)) %>%
      dplyr::mutate_all(list(~na_if(., "Unknown")))
  }

  dfall = data.table::rbindlist(df_list) %>% as.data.frame()
  dfall$nRb = .biohelper_missing_taxonomy_count(dfall, ranks)
  #dfall <- dfall %>% mutate_all(na_if,"")
  dfall <- dfall %>% dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(.x, "")))
  dfall = dfall[!dfall$nRb==length(ranks),]

  if (!nrow(dfall)) {
    empty <- data.frame(feature_id = character(), stringsAsFactors = FALSE)
    for (rank in tolower(ranks)) empty[[rank]] <- character()
    if (ambiguity != "first") attr(empty, "normalisation") <- reports
    return(empty)
  }

  newdf = data.frame(matrix(nrow=dfall$feature_id %>% unique() %>% length(), ncol = length(ranks)+1))
  colnames(newdf) = c("feature_id",ranks)
  newdf$feature_id = dfall$feature_id %>% unique()
  feature_rows <- split(seq_len(nrow(dfall)), factor(dfall$feature_id, levels = newdf$feature_id))

  show_progress <- !use_blast || isTRUE(blast_options$verbose)
  if (show_progress) cat("\nMerging taxonomy for", nrow(newdf), "features\n")
  pb <- if (show_progress) utils::txtProgressBar(min = 0, max = nrow(newdf), initial = 0, style = 3) else NULL
  for (i in 1:nrow(newdf)) {
    # df per ASV
    temp = dfall[feature_rows[[i]], , drop = FALSE]
    # If more than 1 taxo assignment result for this ASV...
    if(temp$df %>% unique() %>% length() >1){
      # If a specific df should be used in case of discrepancy...
      if(!is.na(priority_df)){
        # Assign taxonomy from that df
        temp = temp[temp$df==as.character(priority_df)]
        newdf[i,c("feature_id", ranks)] = temp[which.min(temp$nRb), c("feature_id", ranks), drop = FALSE]
      }else{
        # Determine the number of different values at each rank
        x = vapply(temp[, ranks, drop = FALSE], .biohelper_n_distinct_non_missing, integer(1))
        # While some ranks have more than 1 different values...
        while (!all(x %in% c(0,1))) {
          # Determine the lowest resolution rank where there is discrepancy
          r = min(which(x > 1L))
          r_name <- ranks[[r]]
          # If one of the assignments at that rank has majority...
          rank_table <- table(temp[[r_name]])
          if(length(rank_table) > 0 && max(rank_table) > nrow(temp)/2){
            # Assign the complete taxonomy from the assignment that has majority
            temp = temp[which(temp[[r_name]] == names(which.max(rank_table))),]
            # Determine the number of different values at each rank
            x = vapply(temp[, ranks, drop = FALSE], .biohelper_n_distinct_non_missing, integer(1))
          }else{
            # If no assignment has majority, replace that assignment by NA
            temp[,ranks[r:length(ranks)]] = NA
            # Determine the number of different values at each rank
            x = vapply(temp[, ranks, drop = FALSE], .biohelper_n_distinct_non_missing, integer(1))
          }
        }
      }
      # Once consensus is reached, assigned taxonomy from the assignment with highest resolution
      newdf[i,c("feature_id", ranks)] = temp[which.min(temp$nRb), c("feature_id", ranks), drop = FALSE]
      # In case there is already consensus across all assignment, keep the one with highest resolution
    }else{
      newdf[i,c("feature_id", ranks)] = temp[which.min(temp$nRb), c("feature_id", ranks), drop = FALSE]
    }

    if (show_progress) utils::setTxtProgressBar(pb,i)
  }
  if (show_progress) close(pb)
  # Determine the number of rank assigned for each ASV
  newdf$nRb =  .biohelper_taxonomy_resolution_count(newdf, ranks)
  # Provide the mean and sd of number of rank assigned across the consensus taxo data
  temp_summary = newdf %>% dplyr::summarise(mean = round(mean(nRb),2), sd = round(stats::sd(nRb),2))
  if (show_progress) cat("\nMean assigned taxonomic ranks: ",temp_summary$mean %>% as.numeric(),"\nStandard deviation: ",temp_summary$sd %>% as.numeric(),"\n\n")
  newdf <- newdf %>% dplyr::select(-nRb)
  colnames(newdf) <- tolower(colnames(newdf))
  if (ambiguity != "first") attr(newdf, "normalisation") <- reports
  return(newdf)
}
