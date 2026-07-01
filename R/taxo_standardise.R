.biohelper_taxonomy_output_ranks <- function() {
  c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}

.biohelper_taxonomy_split_ranks <- function(n_fields = 9L) {
  if (n_fields >= 9L) {
    c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  } else {
    c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  }
}

.biohelper_taxonomy_normalisation_ranks <- function(ranks) {
  ranks <- .biohelper_standardise_rank_names(ranks)
  out <- tolower(ranks)
  if ("domain" %in% out) {
    out <- c("domain", "superkingdom", setdiff(out, "domain"))
  }
  unique(out)
}

.biohelper_standardise_rank_names <- function(ranks) {
  ranks <- trimws(as.character(ranks))
  ranks_lower <- tolower(ranks)
  ranks_lower[ranks_lower %in% c("domain", "superkingdom")] <- "domain"
  valid <- tolower(.biohelper_taxonomy_output_ranks())
  ranks_lower <- ranks_lower[ranks_lower %in% valid]
  ranks_lower <- unique(ranks_lower)
  stats::setNames(.biohelper_taxonomy_output_ranks(), tolower(.biohelper_taxonomy_output_ranks()))[ranks_lower]
}

.biohelper_detect_column <- function(df, candidates) {
  column_key <- tolower(colnames(df))
  candidates <- tolower(candidates)
  index <- match(candidates, column_key)
  index <- index[!is.na(index)]
  if (length(index) == 0) {
    return(NA_character_)
  }
  colnames(df)[index[[1]]]
}

.biohelper_clean_taxonomy_value <- function(value) {
  value <- as.character(value)
  value <- trimws(value)
  value[value %in% c("", "NA", "Na", "na", "N/A", "n/a", "Unknown", "unknown")] <- NA_character_
  value
}

.biohelper_coalesce_taxonomy_columns <- function(df, candidates) {
  candidates <- candidates[!is.na(candidates)]
  if (length(candidates) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  out <- rep(NA_character_, nrow(df))
  for (column in candidates) {
    value <- .biohelper_clean_taxonomy_value(df[[column]])
    replace <- is.na(out) & !is.na(value)
    out[replace] <- value[replace]
  }
  out
}

.biohelper_add_split_taxonomy_columns <- function(df, taxonomy_col = "taxonomy") {
  taxonomy_col <- .biohelper_detect_column(df, taxonomy_col)
  if (is.na(taxonomy_col)) {
    return(df)
  }

  taxonomy <- as.character(df[[taxonomy_col]])
  split_taxonomy <- strsplit(taxonomy, ";", fixed = TRUE)
  n_fields <- max(lengths(split_taxonomy), 0L)
  if (n_fields == 0L) {
    return(df)
  }

  split_ranks <- .biohelper_taxonomy_split_ranks(n_fields)
  n_fields <- min(n_fields, length(split_ranks))
  for (rank_index in seq_len(n_fields)) {
    rank <- split_ranks[[rank_index]]
    if (rank %in% tolower(colnames(df))) {
      next
    }
    df[[rank]] <- vapply(
      split_taxonomy,
      function(x) {
        if (length(x) < rank_index) {
          return(NA_character_)
        }
        .biohelper_clean_taxonomy_value(x[[rank_index]])
      },
      character(1)
    )
  }

  df
}

.biohelper_standardise_taxonomy_columns <- function(
  df,
  id_col = NULL,
  output_id_col = NULL,
  taxonomy_col = "taxonomy",
  assignment_method = NULL,
  keep_assignment_method = FALSE
) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df <- .biohelper_add_split_taxonomy_columns(df, taxonomy_col = taxonomy_col)

  if (is.null(id_col)) {
    id_col <- .biohelper_detect_column(
      df,
      c("ASV", "ASVs", "feature_id", "feature.id", "OTU", "OTUs")
    )
  } else {
    id_col <- .biohelper_detect_column(df, id_col)
  }
  if (is.na(id_col)) {
    stop("Could not identify the feature/ASV identifier column.", call. = FALSE)
  }
  if (is.null(output_id_col)) {
    output_id_col <- id_col
  }

  rank_aliases <- list(
    Domain = c("Domain", "domain", "Superkingdom", "superkingdom"),
    Kingdom = c("Kingdom", "kingdom"),
    Phylum = c("Phylum", "phylum"),
    Class = c("Class", "class"),
    Order = c("Order", "order"),
    Family = c("Family", "family"),
    Genus = c("Genus", "genus"),
    Species = c("Species", "species")
  )

  out <- list()
  out[[output_id_col]] <- as.character(df[[id_col]])

  if (isTRUE(keep_assignment_method)) {
    if (!is.null(assignment_method)) {
      out$assignment_method <- rep(as.character(assignment_method), nrow(df))
    } else {
      method_col <- .biohelper_detect_column(df, c("assignment_method", "assignment_source", "method", "source"))
      out$assignment_method <- if (is.na(method_col)) {
        NA_character_
      } else {
        as.character(df[[method_col]])
      }
    }
  }

  for (rank in .biohelper_taxonomy_output_ranks()) {
    columns <- vapply(
      rank_aliases[[rank]],
      function(alias) .biohelper_detect_column(df, alias),
      character(1)
    )
    out[[rank]] <- .biohelper_coalesce_taxonomy_columns(df, columns)
  }

  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

.biohelper_taxonomy_resolution_count <- function(df, ranks = .biohelper_taxonomy_output_ranks()) {
  rowSums(!is.na(df[, ranks, drop = FALSE]) & df[, ranks, drop = FALSE] != "")
}

.biohelper_missing_taxonomy_count <- function(df, ranks = .biohelper_taxonomy_output_ranks()) {
  rowSums(is.na(df[, ranks, drop = FALSE]) | df[, ranks, drop = FALSE] == "")
}

.biohelper_n_distinct_non_missing <- function(value) {
  length(unique(value[!is.na(value) & value != ""]))
}

.biohelper_merge_blast_taxonomy_tables <- function(tables) {
  if (is.null(names(tables)) || any(names(tables) == "")) {
    stop("`tables` must be a named list of BLAST taxonomy tables.", call. = FALSE)
  }

  ranks <- .biohelper_taxonomy_output_ranks()
  blast <- do.call(
    rbind,
    lapply(names(tables), function(method) {
      .biohelper_standardise_taxonomy_columns(
        tables[[method]],
        id_col = "ASVs",
        output_id_col = "ASV",
        assignment_method = method,
        keep_assignment_method = TRUE
      )
    })
  )
  rownames(blast) <- NULL
  blast$nRb <- .biohelper_missing_taxonomy_count(blast, ranks)

  asvs <- unique(blast$ASV)
  newdf <- data.frame(ASV = asvs, stringsAsFactors = FALSE, check.names = FALSE)
  for (rank in ranks) {
    newdf[[rank]] <- NA_character_
  }

  for (i in seq_len(nrow(newdf))) {
    temp_blast <- blast[blast$ASV == newdf$ASV[[i]], , drop = FALSE]
    best <- temp_blast[which.min(temp_blast$nRb), , drop = FALSE]
    rank_counts <- vapply(temp_blast[, ranks, drop = FALSE], .biohelper_n_distinct_non_missing, integer(1))
    conflict <- which(rank_counts > 1L)
    values <- as.character(best[1, ranks, drop = TRUE])
    if (length(conflict) > 0L) {
      values[seq.int(conflict[[1]], length(ranks))] <- NA_character_
    }
    newdf[i, ranks] <- values
  }

  newdf
}

.biohelper_format_blastn_taxo_assignment_output <- function(out) {
  out <- .biohelper_standardise_taxonomy_columns(
    out,
    id_col = NULL,
    output_id_col = "ASV"
  )
  .biohelper_validate_blastn_taxo_assignment_output(out)
  out <- janitor::clean_names(out)
  .biohelper_validate_blastn_taxo_assignment_output(out)
  out
}

.biohelper_validate_blastn_taxo_assignment_output <- function(out) {
  domain_col <- .biohelper_detect_column(out, c("Domain", "domain"))
  bad_domain <- !is.na(domain_col) &&
    any(tolower(trimws(as.character(out[[domain_col]]))) %in% c("blastn", "megablast"), na.rm = TRUE)

  if (bad_domain) {
    stop(
      "Internal blastn_taxo_assignment merge error: Domain contains assignment method values. ",
      "This indicates column shifting before writing the taxonomy output. ",
      "Columns are: ",
      paste(names(out), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(out)
}
