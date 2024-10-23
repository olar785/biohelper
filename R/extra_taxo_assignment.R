#' extra_taxo_assignment
#'
#' @description
#' This function takes in a pyloseq object with taxonomy data and currently add the Protozoa and
#' Archaeplastida (excl. Viridiplantae) groups under the Kingdom rank. This script is particularly useful
#' if taxonomy has been normalised using NCBI curated classification and nomenclature database
#' (using the taxo_normalisation function of this package for example) where protozoans and Archaeplastida (excl. Viridiplantae)
#' are generally not assigned any information under Kingdom, making it difficult to dissociate metazoans from protists, or the Plantae group.
#' Importantly, this function requires the download of NCBI taxonomic database.
#' This takes few minutes to install using the following command:\cr
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package. The database is ~65 GB but the user can set getAccessions=FALSE to drastically reduce its size.
#'
#' @param
#' physeq_obj         Phyloseq object
#' @param
#' sqlFile_path       Path to the local NCBI taxonomy db
#'
#'
#' @export
#' @examples
#' data("ps_test_data")
#' extra_taxo_assignment(ps_test_data, sqlFile = 'accessionTaxa.sql')

extra_taxo_assignment = function(physeq_obj, sqlFile_path){
  # Taxonomic ID for protists
  proto_names <- c("Telonemia", "Stramenopiles", "Alveolata", "Rhizaria","Chromeraceae", # TSAR
                   "Tubulinea", "Discosea", "Evosea","Elardia","Flabellinia","Echinamoebida","Breviatea","Filasterea", # Amoebozoa
                   "Euglenozoa", "Kinetoplastea", # Discoba
                   "Hemimastigophora","Nibbleridia","Centroplasthelida","Ichthyosporea","Nibbleridea","Ancyromonadida","Jakobida", # ?
                   "Heterolobosea", # Diphoda
                   "Choanoflagellata") # Choanozoa
  other_proto = c("Telonemida","Chromeraceae","Apusomonadidae","Apusomonadida","Nucleariidae","Rotosphaerida","Colpodellida","Colpodellaceae") # May not be in NCBI db
  archaeplastida_names <- c("Bangiophyceae","Chrysophyceae","Florideophyceae","Cryptophyceae","Rhodophyta","Haptophyta","Picozoa","Picomonadea")
  aphelidea_names <- c("Aphelidea") #fungi

  # Step 1: Use getId() to get the taxonomic ID for each taxon name
  proto_ids <- taxonomizr::getId(taxa = proto_names, sqlFile = sqlFile_path)
  archaeplastida_ids <- taxonomizr::getId(taxa = archaeplastida_names, sqlFile = sqlFile_path)
  aphelidea_ids <- taxonomizr::getId(taxa = aphelidea_names, sqlFile = sqlFile_path)

  # Step 2: Get all descendant taxa for TSAR taxa
  proto_taxa <- taxonomizr::getDescendants(ids = proto_ids, sqlFile = sqlFile_path, desiredTaxa = c("phylum", "class","order"))
  proto_taxa = c(proto_taxa,proto_names)
  archaeplastida_taxa <- taxonomizr::getDescendants(ids = archaeplastida_ids, sqlFile = sqlFile_path, desiredTaxa = c("phylum", "class","order"))
  archaeplastida_taxa = c(archaeplastida_taxa,archaeplastida_names)
  aphelidea_taxa <- taxonomizr::getDescendants(ids = aphelidea_ids, sqlFile = sqlFile_path, desiredTaxa = c("phylum", "class","order"))
  aphelidea_taxa = c(aphelidea_taxa,aphelidea_names)

  # Step 3: Extract the taxonomic table
  tax_table_df <- as.data.frame(phyloseq::tax_table(physeq_obj))
  # Step 4: Define a function to check and update the kingdom level
  update_kingdom <- function(tax_row) {
    # Check if Kingdom is NA or an empty string, and if any of the higher ranks belong to TSAR
    if ((is.na(tax_row["kingdom"]) || tax_row["kingdom"] == "") &&
        (tax_row["phylum"] %in% proto_taxa ||
         tax_row["class"] %in% proto_taxa ||
         tax_row["order"] %in% proto_taxa)) {
      tax_row["kingdom"] <- "Protozoa"
    }
    if ((is.na(tax_row["kingdom"]) || tax_row["kingdom"] == "") &&
        (tax_row["phylum"] %in% other_proto ||
         tax_row["class"] %in% other_proto ||
         tax_row["order"] %in% other_proto ||
         tax_row["family"] %in% other_proto)) {
      tax_row["kingdom"] <- "Protozoa"
    }
   if ((is.na(tax_row["kingdom"]) || tax_row["kingdom"] == "") &&
           (tax_row["phylum"] %in% archaeplastida_taxa ||
            tax_row["class"] %in% archaeplastida_taxa ||
            tax_row["order"] %in% archaeplastida_taxa)) {
    tax_row["kingdom"] <- "Archaeplastida (excl. Viridiplantae)"
   }
    if ((is.na(tax_row["kingdom"]) || tax_row["kingdom"] == "") &&
        (tax_row["phylum"] %in% aphelidea_taxa ||
         tax_row["class"] %in% aphelidea_taxa ||
         tax_row["order"] %in% aphelidea_taxa)) {
      tax_row["kingdom"] <- "Fungi"
    }
    return(tax_row)
  }

  # Step 5: Apply the update function to all rows in the tax_table_df
  tax_table_df <- as.data.frame(t(apply(tax_table_df, 1, update_kingdom)))

  # Step 6: Reassign the modified tax_table to the phyloseq object
  phyloseq::tax_table(physeq_obj) <- phyloseq::tax_table(as.matrix(tax_table_df))
  return(physeq_obj)
}

