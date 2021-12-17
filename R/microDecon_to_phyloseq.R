#' Load a micronDecon_2_phyloseq
#'
#' This function performs normalisation of taxonomic assignments of 18S and COI data
#' identified via silva, midori or any other reference databases, using NCBI curated
#' classification and nomenclatureand database. It can also provide further taxonomic
#' information of additional desired ranks.
#'
#' The first column must contain ASV/OTU identifiers as the example below:
#' ASV     Kingdom   Phylum   etc
#' ASV_1   Eukaryota   Arthropoda
#' ASV_2   Eukaryota   Chordata
#' ASV_3   Eukaryota   Chordata
#'
#' Importantly, this function requires the download of NCBI taxonomic database (~65 GB).
#' This takes few minutes to install using the following command:
#' prepareDatabase('accessionTaxa.sql') # From the taxonomizr R package

#'
#' @export
#' @examples
#' micronDecon_2_phyloseq(taxo_dataframe, sqlFile = sqlFile, desired_ranks = c("Family", "Genus", "Species"))

# 1 - Decontamination of DNA extraction blanks per batch (extraction type)
# 2 - Decontamination of PCR and indexing/sequencing blanks overall
microDecon_to_phyloseq = function(ps, env, decontaminated){
  ps = prune_samples(sample_sums(ps) > 0, ps)
  ps@sam_data$amplicon_type = ps@sam_data$amplicon_type %>% tolower()
  # 1.1 - Decontamination of DNA extraction blanks per batch (extraction type).
  ps_ext = ps %>% subset_samples(amplicon_type != "pcr blank" & amplicon_type != "indexing/sequencing blank") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  metadata_ext = pstoveg_sample(ps_ext) %>% dplyr::arrange(amplicon_type)
  ASV_table_ps = pstoveg_otu(ps_ext) %>% t() %>% as.data.frame()
  ASV_table_ps = ASV_table_ps[,match(metadata_ext$sample_name, names(ASV_table_ps))] %>% rownames_to_column('OTU_ID')

  names_blanks_ps = c("OTU_ID",ps_ext %>% subset_samples(type != "sample") %>% sample_names())
  names_samples_ps = ps_ext %>% subset_samples(type == "sample") %>% sample_names()
  ASV_table_ps <- subset(ASV_table_ps, select=c(names_blanks_ps,names_samples_ps))

  Taxo=as.data.frame(ps_ext@tax_table@.Data)
  Taxo$OTU_ID = rownames(Taxo)
  sorted <- Taxo %>% arrange(Kingdom,Phylum,Class,Order,Family,Genus,Species)
  sorted2 = sorted %>% unite("Taxonomy", Kingdom:Species,sep = ";")
  ASV_table_ps = ASV_table_ps[match(sorted2$OTU_ID,ASV_table_ps$OTU_ID),]
  #all(ASV_table_ps$OTU_ID==sorted2$OTU_ID)
  ASV_table_ps$Taxonomy = sorted2$Taxonomy

  df_temp = metadata_ext %>% group_by(type) %>% dplyr::count()
  #MicroDecon function
  decontaminated_ext_1 <- decon(data = ASV_table_ps, numb.blanks=df_temp$n[1], numb.ind= c(df_temp$n[2]),taxa = TRUE)
  ps_trimmed_1.1 = micronDecon_2_phyloseq(ps_18S = ps, env = pstoveg_sample(ps), decontaminated_ext_1)
  ps_out = ps_18S %>% subset_samples(sample_names(ps_18S) %ni% sample_names(ps_trimmed_1.1)) %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  ps_trimmed_c = merge_phyloseq(ps_out, ps_trimmed_1.1)






  if("Mean.blank" %in% colnames(decontaminated$decon.table)){
    otu_table_ps = otu_table(decontaminated$decon.table %>% dplyr::select(-c(OTU_ID,Mean.blank,Taxonomy)), taxa_are_rows = T)
  }else{
    otu_table_ps = otu_table(decontaminated$decon.table %>% dplyr::select(-c(OTU_ID,Taxonomy)), taxa_are_rows = T)
  }
  rownames(otu_table_ps) = decontaminated$decon.table$OTU_ID
  #tax_ps = tax_table(decontaminated$decon.table$Taxonomy %>% colsplit(";", names = c("Kingdom","phylum","class","order","family","genus","ASV")))
  tax_ps = tax_table(decontaminated$decon.table$Taxonomy %>% colsplit(";", names = c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")))
  rownames(tax_ps) = decontaminated$decon.table$OTU_ID
  ps_trimmed = merge_phyloseq(otu_table_ps,tax_ps)
  # Adding back the sequences
  ps_taxa_trimmed = prune_taxa(ps_trimmed %>% taxa_names(), ps_18S)
  fasta_ASVs = ps_taxa_trimmed@refseq
  fasta_ASVs = fasta_ASVs[match(ps_taxa_trimmed %>% taxa_names(), fasta_ASVs@ranges@NAMES),]
  #all(fasta_ASVs@ranges@NAMES == ps_trimmed %>% taxa_names())
  ps_trimmed@refseq = phyloseq::refseq(fasta_ASVs)
  # Adding back the environmental data
  env=env[rownames(env) %in% sample_names(ps_trimmed),]
  env=env[match(sample_names(ps_trimmed),rownames(env)),]
  #all(env$ngs_id==sample_names(ps_trimmed))
  sample_data(ps_trimmed) = sample_data(env)
  ps_trimmed
  #colnames(ps_trimmed@tax_table) = c("Kingdom","phylum", "class","order","family","genus","species","ASV")
  colnames(ps_trimmed@tax_table) = c("Kingdom","Phylum", "Class","Order","Family","Genus","Species")
  return(ps_trimmed)
}
