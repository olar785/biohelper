.biohelper_managed_ncbi_path <- function() {
  directory <- if (exists("R_user_dir", asNamespace("tools"), inherits = FALSE)) {
    tools::R_user_dir("biohelper", "cache")
  } else file.path(path.expand("~"), ".cache", "biohelper")
  file.path(directory, "ncbi.sqlite")
}

.biohelper_lca_db_path <- function(taxonomy_db = NULL) {
  if (!is.null(taxonomy_db)) {
    if (!is.character(taxonomy_db) || length(taxonomy_db) != 1L || is.na(taxonomy_db) || !nzchar(taxonomy_db)) stop("`taxonomy_db` must be a single SQLite file path.", call. = FALSE)
    return(path.expand(taxonomy_db))
  }
  managed <- .biohelper_managed_ncbi_path()
  if (file.exists(managed)) managed else path.expand(taxizedb::db_path("ncbi"))
}

.biohelper_ncbi_schema <- function(connection) {
  tables <- DBI::dbListTables(connection)
  if (!all(c("names", "nodes") %in% tables)) stop("Required NCBI tables 'names' and 'nodes' are absent.")
  node_columns <- DBI::dbListFields(connection, "nodes")
  name_columns <- DBI::dbListFields(connection, "names")
  if (all(c("id", "parent", "rank") %in% node_columns) && all(c("id", "name", "scientific") %in% name_columns)) {
    list(node_id = "id", parent = "parent", name_id = "id", name = "name", scientific = "scientific = 1")
  } else if (all(c("tax_id", "parent_tax_id", "rank") %in% node_columns) && all(c("tax_id", "name_txt", "name_class") %in% name_columns)) {
    list(node_id = "tax_id", parent = "parent_tax_id", name_id = "tax_id", name = "name_txt", scientific = "name_class = 'scientific name'")
  } else stop("Unsupported or incomplete NCBI names/nodes schema.")
}

.biohelper_validate_ncbi_db <- function(path, thorough = FALSE) {
  path <- path.expand(path)
  if (!file.exists(path) || dir.exists(path)) stop("NCBI taxonomy database is missing: ", path, ". Supply `taxonomy_db` or run with update=TRUE to build a separate biohelper cache.", call. = FALSE)
  if (file.access(path, 4L) != 0L) stop("NCBI taxonomy database is inaccessible (read permission): ", path, call. = FALSE)
  connection <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO, synchronous = NULL),
    error = function(e) stop("Cannot open NCBI taxonomy database: ", path, ". ", conditionMessage(e), call. = FALSE))
  on.exit(DBI::dbDisconnect(connection))
  tryCatch({
    schema <- .biohelper_ncbi_schema(connection)
    probe <- DBI::dbGetQuery(connection, paste0(
      "SELECT n.", schema$node_id, " AS id, n.", schema$parent, " AS parent, n.rank, s.", schema$name,
      " AS name FROM nodes n JOIN names s ON n.", schema$node_id, " = s.", schema$name_id,
      " WHERE s.", schema$scientific, " LIMIT 1"))
    if (nrow(probe) != 1L || anyNA(probe[, c("id", "parent", "name")])) stop("No usable scientific-name/node records.")
    if (thorough) {
      check <- DBI::dbGetQuery(connection, "PRAGMA quick_check")[[1L]]
      if (!identical(check, "ok")) stop(paste(check, collapse = "; "))
      broken <- DBI::dbGetQuery(connection, paste0("SELECT n.", schema$node_id,
        " FROM nodes n LEFT JOIN nodes p ON n.", schema$parent, " = p.", schema$node_id,
        " WHERE n.", schema$parent, " != 1 AND p.", schema$node_id, " IS NULL LIMIT 1"))
      if (nrow(broken)) stop("Downloaded taxonomy contains missing parent nodes.")
      unnamed <- DBI::dbGetQuery(connection, paste0("SELECT n.", schema$node_id,
        " FROM nodes n LEFT JOIN names s ON n.", schema$node_id, " = s.", schema$name_id,
        " AND s.", schema$scientific, " WHERE s.", schema$name_id, " IS NULL LIMIT 1"))
      if (nrow(unnamed)) stop("Downloaded taxonomy contains nodes without scientific names.")
    }
    list(path = normalizePath(path), schema = schema,
         identity = .biohelper_file_fingerprint(path, content = FALSE),
         provenance = if ("biohelper_provenance" %in% DBI::dbListTables(connection)) DBI::dbReadTable(connection, "biohelper_provenance") else NULL)
  }, error = function(e) stop("NCBI taxonomy database is corrupt or not functional: ", path, ". ", conditionMessage(e), call. = FALSE))
}

.biohelper_download_ncbi_archive <- function(destination) {
  url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
  if (requireNamespace("curl", quietly = TRUE)) curl::curl_download(url, destination, quiet = FALSE)
  else utils::download.file(url, destination, mode = "wb", method = "libcurl", quiet = FALSE)
  invisible(url)
}

.biohelper_build_ncbi_archive <- function(archive, output, workdir) {
  magic <- readBin(archive, "raw", n = 4L)
  if (length(magic) != 4L || !identical(magic[1:2], charToRaw("PK"))) stop("Downloaded NCBI archive is not a valid ZIP (unexpected file header).", call. = FALSE)
  members <- tryCatch(suppressWarnings(utils::unzip(archive, list = TRUE)), error = function(e) stop("Downloaded NCBI archive is not a valid ZIP: ", conditionMessage(e), call. = FALSE))
  required <- c("names.dmp", "nodes.dmp")
  if (!all(required %in% members$Name)) stop("Downloaded NCBI ZIP does not contain names.dmp and nodes.dmp.", call. = FALSE)
  wanted <- intersect(c(required, "merged.dmp"), members$Name)
  withCallingHandlers(utils::unzip(archive, files = wanted, exdir = workdir),
    warning = function(w) stop("NCBI ZIP extraction failed: ", conditionMessage(w), call. = FALSE))
  for (name in required) if (!file.exists(file.path(workdir, name)) || file.info(file.path(workdir, name))$size == 0) stop("Incomplete NCBI archive member: ", name, call. = FALSE)
  taxonomizr::read.names.sql(file.path(workdir, "names.dmp"), output)
  taxonomizr::read.nodes.sql(file.path(workdir, "nodes.dmp"), output)
  connection <- DBI::dbConnect(RSQLite::SQLite(), output)
  on.exit(DBI::dbDisconnect(connection))
  if ("merged.dmp" %in% wanted) {
    merged <- data.table::fread(file.path(workdir, "merged.dmp"), sep = "|", header = FALSE, select = 1:2,
                               col.names = c("old_id", "new_id"), showProgress = FALSE)
    DBI::dbWriteTable(connection, "merged", as.data.frame(merged))
    DBI::dbExecute(connection, "CREATE INDEX merged_old_id ON merged(old_id)")
  }
  provenance <- data.frame(format = "biohelper_ncbi_compact_v1",
    source = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",
    built_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    archive_md5 = unname(tools::md5sum(archive)),
    taxonomizr_version = as.character(utils::packageVersion("taxonomizr")), stringsAsFactors = FALSE)
  DBI::dbWriteTable(connection, "biohelper_provenance", provenance)
}

.biohelper_update_ncbi_safe <- function(path) {
  # Never replace a shared taxizedb/taxonomizr database with our compact schema.
  # Existing biohelper-owned databases can be replaced, with an archived copy.
  if (file.exists(path)) {
    existing <- tryCatch(.biohelper_validate_ncbi_db(path), error = function(e) {
      if (identical(normalizePath(path, mustWork = FALSE), normalizePath(.biohelper_managed_ncbi_path(), mustWork = FALSE))) return(NULL)
      stop(e)
    })
    is_managed_path <- identical(normalizePath(path, mustWork = FALSE), normalizePath(.biohelper_managed_ncbi_path(), mustWork = FALSE))
    if (!is_managed_path && (is.null(existing$provenance) || !"format" %in% names(existing$provenance) ||
        !all(existing$provenance$format == "biohelper_ncbi_compact_v1"))) {
      stop("Refusing to replace an externally managed taxonomy database. Use update=TRUE with taxonomy_db=NULL (a separate biohelper cache), or supply a new file path. Existing database preserved: ", path, call. = FALSE)
    }
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  workdir <- tempfile("ncbi_update_", tmpdir = dirname(path))
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE))
  archive <- file.path(workdir, "taxdmp.zip")
  replacement <- file.path(workdir, "ncbi.sqlite")
  message("NCBI taxonomy update: downloading over HTTPS into a separate staging directory")
  .biohelper_download_ncbi_archive(archive)
  .biohelper_build_ncbi_archive(archive, replacement, workdir)
  .biohelper_validate_ncbi_db(replacement, thorough = TRUE)
  .biohelper_atomic_publish(replacement, path, overwrite = FALSE,
                            archive_dir = file.path(dirname(path), "archive"))
  .biohelper_validate_ncbi_db(path)
}

.biohelper_prepare_lca_db <- function(update = FALSE, taxonomy_db = NULL) {
  if (!is.logical(update) || length(update) != 1L || is.na(update)) stop("`update` must be TRUE or FALSE.", call. = FALSE)
  if (update) {
    target <- if (is.null(taxonomy_db)) .biohelper_managed_ncbi_path() else path.expand(taxonomy_db)
    .biohelper_update_ncbi_safe(target)
  } else .biohelper_validate_ncbi_db(.biohelper_lca_db_path(taxonomy_db))
}

.biohelper_lca_lookup <- function(ids, database, verbose = FALSE) {
  ids <- sort(unique(as.character(ids)))
  requested <- ids[!is.na(ids) & grepl("^[1-9][0-9]*$", ids)]
  ranks <- c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  out <- matrix(NA_character_, length(ids), length(ranks), dimnames = list(ids, ranks))
  if (!length(requested)) return(list(taxonomy = out, resolved = 0L, unresolved = length(ids)))
  if (!identical(.biohelper_file_fingerprint(database$path, content = FALSE), database$identity)) stop("Taxonomy database changed after validation; retry with one stable database version.", call. = FALSE)
  connection <- DBI::dbConnect(RSQLite::SQLite(), database$path, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(connection))
  schema <- database$schema
  records <- new.env(parent = emptyenv())
  assign("1", list(parent = "1", rank = "no rank", name = "root"), records)
  pending <- requested[requested != "1"]
  repeat {
    pending <- unique(pending[!vapply(pending, exists, logical(1), envir = records, inherits = FALSE)])
    if (!length(pending)) break
    next_ids <- character()
    for (batch in split(pending, ceiling(seq_along(pending) / 800L))) {
      query <- paste0("SELECT n.", schema$node_id, " AS id, n.", schema$parent,
        " AS parent, n.rank, s.", schema$name, " AS name FROM nodes n LEFT JOIN names s ON n.",
        schema$node_id, " = s.", schema$name_id, " AND s.", schema$scientific,
        " WHERE n.", schema$node_id, " IN (", paste(rep("?", length(batch)), collapse = ","), ")")
      found <- tryCatch(DBI::dbGetQuery(connection, query, params = as.list(as.numeric(batch))),
        error = function(e) stop("Systemic NCBI taxonomy lookup failure: ", conditionMessage(e), call. = FALSE))
      if (anyDuplicated(found$id)) stop("NCBI database contains duplicate scientific records for a TaxID.", call. = FALSE)
      for (id in batch) {
        index <- match(id, as.character(found$id))
        value <- if (is.na(index)) NULL else list(parent = as.character(found$parent[[index]]), rank = found$rank[[index]], name = found$name[[index]])
        assign(id, value, records)
        if (!is.null(value)) {
          if (is.na(value$parent) || is.na(value$name) || is.na(value$rank)) stop("Corrupt NCBI taxonomy: node has no parent, rank or scientific name (TaxID ", id, ").", call. = FALSE)
          next_ids <- c(next_ids, value$parent)
        }
      }
    }
    pending <- next_ids
  }
  resolved <- 0L
  if (!identical(.biohelper_file_fingerprint(database$path, content = FALSE), database$identity)) stop("Taxonomy database changed during lookup; retry with one stable database version.", call. = FALSE)
  for (id in requested) {
    current <- id
    seen <- character()
    first <- TRUE
    while (current != "1") {
      if (current %in% seen || length(seen) >= 200L) stop("Corrupt NCBI taxonomy: cyclic or excessive ancestry at TaxID ", id, call. = FALSE)
      seen <- c(seen, current)
      value <- get(current, records, inherits = FALSE)
      if (is.null(value)) {
        if (!first) stop("Corrupt NCBI taxonomy: missing ancestor ", current, " for TaxID ", id, call. = FALSE)
        break
      }
      rank <- if (value$rank == "superkingdom") "domain" else value$rank
      if (rank %in% ranks && is.na(out[id, rank])) out[id, rank] <- value$name
      current <- value$parent
      first <- FALSE
    }
    if (!first) resolved <- resolved + 1L
  }
  if (verbose) message("NCBI TaxIDs: ", resolved, " resolved; ", length(ids) - resolved, " unresolved")
  list(taxonomy = out, resolved = resolved, unresolved = length(ids) - resolved)
}
