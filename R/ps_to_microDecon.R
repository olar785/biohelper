#' Loads ps_to_microDecon
#'
#' This function is a wrapper of the decon function from microDecon R package.
#' It performs the decon function on a phyloseq object with sample data and
#' returns a decontaminated phyloseq object with sample data, taxonomy and reference sequences if present.
#' The sample_data MUST have a column labeled sample_id and a column labeled amplicon_type
#' Non blank samples must be labeled 'sample' or 'Sample'.
#'
#'For example:
#'sample_id   amplicon_type         DNA_extraction_batch      extraction_method   etc.
#'sample1     sample                1                         manual
#'sample2     sampling_blank        1                         manual
#'sample3     dna_extraction_blank  2                         robot
#'sample4     sample                2                         robot
#'sample5     pcr_blank             2                         robot
#'
#' @export
#' @examples
#' ps_to_microDecon(ps_test_data, groups = "extraction_method")

ps_to_microDecon = function(ps, groups=NA, runs=2, thresh = 0.7, prop.thresh = 0.00005, regression = 0, low.threshold=40, up.threshold=400){
  microDecon_2_phyloseq = function(ps_obj, env, decontaminated, taxo_ranks=NULL){
    if("Mean.blank" %in% colnames(decontaminated$decon.table)){
      otu_table_ps = otu_table(decontaminated$decon.table[colnames(decontaminated$decon.table) %ni% c("OTU_ID","Mean.blank","Taxonomy")], taxa_are_rows = T)
    }else{
      otu_table_ps = otu_table(decontaminated$decon.table[colnames(decontaminated$decon.table) %ni% c("OTU_ID","Taxonomy")], taxa_are_rows = T)
    }
    rownames(otu_table_ps) = decontaminated$decon.table$OTU_ID
    if(!is.null(taxo_ranks)){
      tax_ps = tax_table(decontaminated$decon.table$Taxonomy %>% colsplit(";", names = taxo_ranks) %>% as.matrix())
      rownames(tax_ps) = decontaminated$decon.table$OTU_ID
      ps_trimmed = merge_phyloseq(otu_table_ps,tax_ps)
      colnames(ps_trimmed@tax_table) = taxo_ranks
    }else{
      ps_trimmed = otu_table_ps
    }
    if(!is.null(ps_obj@refseq)){
      # Adding back the sequences
      ps_taxa_trimmed = prune_taxa(ps_trimmed %>% taxa_names(), ps_obj)
      fasta_ASVs = ps_taxa_trimmed@refseq
      fasta_ASVs = fasta_ASVs[match(ps_taxa_trimmed %>% taxa_names(), fasta_ASVs@ranges@NAMES),]
      ps_trimmed = merge_phyloseq(ps_trimmed,phyloseq::refseq(fasta_ASVs))
    }
    # Adding back the environmental data
    env=env[rownames(env) %in% sample_names(ps_trimmed),]
    env=env[match(sample_names(ps_trimmed),rownames(env)),]
    sample_data(ps_trimmed) = sample_data(env)
    return(ps_trimmed %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE))
  }

  ps = prune_samples(sample_sums(ps) > 0, ps) %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  #ps@sam_data$amplicon_type = suppressWarnings(str_replace_all(ps@sam_data$amplicon_type, pattern = c("NA","na","Na","NaN","nan",""), replacement = NA_character_))
  ps@sam_data$amplicon_type = ps@sam_data$amplicon_type %>% tolower()
  metadata = pstoveg_sample(ps) %>% dplyr::arrange(amplicon_type)
  ASV_table_ps = pstoveg_otu(ps) %>% t() %>% as.data.frame()
  ASV_table_ps = ASV_table_ps[,match(metadata$sample_id, names(ASV_table_ps))] %>% rownames_to_column('OTU_ID')
  names_blanks_ps = c("OTU_ID",ps %>% subset_samples(amplicon_type!="sample") %>% sample_names())
  names_samples_ps = ps %>% subset_samples(amplicon_type=="sample") %>% sample_names()
  ASV_table_ps <- subset(ASV_table_ps, select=c(names_blanks_ps,names_samples_ps))
  if(!is.null(ps@tax_table)){
    Taxo=as.data.frame(ps@tax_table@.Data); Taxo$OTU_ID = rownames(Taxo)
    sorted = Taxo %>% unite("Taxonomy", 1:(ncol(Taxo)-1),sep = ";")
    ASV_table_ps = ASV_table_ps[match(sorted$OTU_ID,ASV_table_ps$OTU_ID),]
    ASV_table_ps$Taxonomy = sorted$Taxonomy
  }
  df_temp = metadata %>% group_by_at(which(colnames(metadata) %in% c("amplicon_type",groups))) %>% dplyr::count()
  #MicroDecon function
  if(!is.null(ps@tax_table)){
    decontaminated_ext <- decon(data = ASV_table_ps, numb.blanks=sum(df_temp[df_temp$amplicon_type!="sample",]$n), numb.ind = df_temp[df_temp$amplicon_type=="sample",]$n, taxa = TRUE,runs, thresh, prop.thresh, regression, low.threshold, up.threshold)
    ps_trimmed = microDecon_2_phyloseq(ps_obj = ps, env = pstoveg_sample(ps), decontaminated = decontaminated_ext, taxo_ranks = colnames(Taxo)[(colnames(Taxo)!="OTU_ID")])
  }else{
    decontaminated_ext <- decon(data = ASV_table_ps, numb.blanks=sum(df_temp[df_temp$amplicon_type!="sample",]$n), numb.ind = df_temp[df_temp$amplicon_type=="sample",]$n, taxa = FALSE,runs, thresh, prop.thresh, regression, low.threshold, up.threshold)
    ps_trimmed = microDecon_2_phyloseq(ps_obj = ps, env = pstoveg_sample(ps), decontaminated = decontaminated_ext)
  }
  # Printing results
  ntaxa_before = ps %>% subset_samples(amplicon_type=="sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% ntaxa()
  ntaxa_after = ntaxa(ps_trimmed)

  cat("\nContamination removal outcome with MicroDecon\n")
  cat(paste0("Number of ASVs totally removed: ",ntaxa_before - ntaxa_after))
  cat("\n")
  cat(paste0("Percent of ASVs removed: ",round((1 - (ntaxa_after / ntaxa_before)) * 100,2), " %"))
  cat("\n")

  reads_before = ps %>% subset_samples(amplicon_type=="sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% sample_sums() %>% sum()
  reads_after = ps_trimmed %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% sample_sums() %>% sum()

  cat(paste0("Total number of reads  removed: ",reads_before - reads_after))
  cat("\nPercent of reads removed: ",paste0(round((1 - (reads_after / reads_before)) * 100,2)," %"))
  cat("\n")
  return(ps_trimmed)
}
#
