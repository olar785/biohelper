# Small synthetic taxonomy for controlled hierarchy/threshold tests. This is
# separate from the authentic NCBI fixture used by integration tests.
workflow_database <- function(schema = "compact") {
  path <- tempfile(fileext = ".sqlite")
  connection <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(connection))
  nodes <- data.frame(id = 1:15,
    parent = c(1, 1, 2, 3, 4, 5, 6, 7, 8, 7, 10, 2, 12, 13, 14),
    rank = c("no rank", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species",
             "genus", "species", "kingdom", "phylum", "genus", "species"))
  names <- data.frame(id = 1:15,
    name = c("root", "Eukaryota", "Metazoa", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "Homo sapiens",
             "Pan", "Pan troglodytes", "Viridiplantae", "Streptophyta", "PlantGenus", "Plant species"), scientific = 1L)
  if (schema == "taxizedb") {
    names(nodes)[1:2] <- c("tax_id", "parent_tax_id")
    names(names) <- c("tax_id", "name_txt", "scientific")
    names$name_class <- "scientific name"
    names$scientific <- NULL
  }
  DBI::dbWriteTable(connection, "nodes", nodes)
  DBI::dbWriteTable(connection, "names", names)
  path
}

workflow_hits <- function(ids = "asv1", taxids = "9", identity = 99, evalue = 1e-30,
                          bitscore = 200, coverage = 100, subject = "reference") {
  data.frame(qseqid = ids, qlen = 100, pident = identity, sseqid = subject, sgi = "0",
             evalue = evalue, bitscore = bitscore, staxids = as.character(taxids),
             sscinames = "Reference organism", sblastnames = "organism", qcovs = coverage,
             qcovhsp = coverage, stringsAsFactors = FALSE)
}

workflow_write_hits <- function(hits, path = tempfile(fileext = ".csv")) {
  if (nrow(hits)) utils::write.table(hits[, biohelper:::.biohelper_blast_raw_columns()], path,
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) else file.create(path)
  path
}

workflow_setup <- function(ids = c("asv1", "asv2")) {
  directory <- tempfile("workflow path ")
  dir.create(directory)
  queries <- file.path(directory, "query file.fasta")
  writeLines(as.vector(rbind(paste0(">", ids), rep("ACGTACGT", length(ids)))), queries)
  database <- workflow_database()
  list(queries = queries, taxonomy_db = database,
       output_path = file.path(directory, "results path"), blastapp_path = "/mock path/blastn",
       db = "/mock database/nt", nthreads = 1, method = "both")
}

workflow_taxonomy <- function(result) {
  biohelper:::.biohelper_standardise_taxonomy_columns(result, id_col = "ASVs", output_id_col = "ASV")
}
