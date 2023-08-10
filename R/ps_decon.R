#' ps_decon
#'
#' @description
#' This function offers 3 different options to remove putative sequence contamination
#' of a phyloseq object. The microDecon option is a wrapper of the decon function
#' from microDecon R package. It performs the decon function on a phyloseq object
#' with sample data and returns a decontaminated phyloseq object with sample data,
#' taxonomy and reference sequences if present.
#'
#' The 'max_v' option subtract read associated to putative contaminant ASV by using
#' their max read count in blank(s).
#'
#' Finally, the 'complete_asv_removal' removes all ASVs found all blanks from the dataset.
#'
#' Note that the 'max_v' option does not take into account the compositional nature of the data.
#'
#' Also note that microDecon assumes a common source of contamination. Quoting the authors,
#' "If substantial differences among blanks occur only across experimental blocks, such as
#' extraction kits (suggesting consistent contamination within a block), then use microDecon
#' separately for each block. If, however, there is substantial variability among blanks within
#' blocks (suggesting contamination from poor laboratory techniques), microDecon will not be effective.
#'
#' The sample_data must have a column labeled sample_id and a column labeled amplicon_type, and the
#' non blank samples must be labeled as 'sample'. Alternatively, the user can provide a list of samples
#' to use as controls.
#'
#' If you use the decon function for a publication, please cite the creators of the microDecon R package (citation('microDecon')).

#'
#'For example:\cr
#' \tabular{rrrrr}{
#'   \strong{sample_id} \tab \strong{amplicon_type} \tab \strong{extraction_batch} \tab \strong{extraction_method} \tab \strong{etc.} \cr
#'   sample1 \tab sample \tab 1 \tab manual \tab NA\cr
#'   sample2 \tab sampling_blank \tab 1 \tab manual \tab NA\cr
#'   sample3 \tab dna_extraction_blank \tab 2 \tab  robot \tab NA\cr
#'   sample4 \tab sample \tab 2 \tab robot \tab NA\cr
#'   sample5 \tab pcr_blank \tab 2 \tab robot \tab NA
#' }
#'
#'
#' @param
#' ps Phyloseq object to decontaminate\cr
#' @param
#' method Method to be used for decontamination. Options are 'microDecon' (using the decon function of microDecon), 'max_v' and 'complete_asv_removal'\cr
#' @param
#' group Will be used in the numb.ind argument of the microDecon::decon function\cr
#' @param
#' (...) If using microDecon the user can specify any argument of the decon function with the exception of num.blanks and numb.ind, which are already handled by ps_decon.\cr
#'
#' @export
#' @examples
#' ps_decon(ps_test_data, method = "microDecon", group = "extraction_method")

ps_decon = function(ps, method= "complete_asv_removal", group=NA, runs=2, thresh = 0.7, prop.thresh = 0.00005, regression = 0, low.threshold=40, up.threshold=400){
  # Creating custom microDecon function
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
    ps_trimmed = ps_trimmed %>%
      phyloseq::subset_samples(amplicon_type=="sample") %>%
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
    return(ps_trimmed)
  }

  # Ensuring sample_id and amplicon_type are present in the metadata
  if(is.null(ps@sam_data$amplicon_type)){
    cat("\n")
    cat("The 'amplicon_type' column is missing from the metadata. Please indicate which rows are 'samples' and which are 'blanks' under 'amplicon_type'.","\n")
    stop()
  }
  if(is.null(ps@sam_data$sample_id)){
    cat("\n")
    cat("The 'sample_id' column is missing from the metadata. Adding it to the phyloseq object using sample_names.","\n")
    ps@sam_data$sample_id = sample_names(ps)
    stop()
  }

  # Ensuring no empty sample exist
  ps = prune_samples(sample_sums(ps) > 0, ps) %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  ASVs_in_Blanks = ps %>% phyloseq::subset_samples(amplicon_type != "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% phyloseq::taxa_names()

  # Processing according to the chosen method
  if(method == "microDecon"){
    #ps@sam_data$amplicon_type = suppressWarnings(str_replace_all(ps@sam_data$amplicon_type, pattern = c("NA","na","Na","NaN","nan",""), replacement = NA_character_))
    ps@sam_data$amplicon_type = ps@sam_data$amplicon_type %>% tolower()
    metadata = theseus::pstoveg_sample(ps) %>% dplyr::arrange(amplicon_type)
    ASV_table_ps = theseus::pstoveg_otu(ps) %>% t() %>% as.data.frame()
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
    df_temp = metadata %>% group_by_at(which(colnames(metadata) %in% c("amplicon_type",group))) %>% dplyr::count()
    #MicroDecon function
    if(!is.null(ps@tax_table)){
      decontaminated_ext <- microDecon::decon(data = ASV_table_ps, numb.blanks=sum(df_temp[df_temp$amplicon_type!="sample",]$n), numb.ind = df_temp[df_temp$amplicon_type=="sample",]$n, taxa = TRUE,runs, thresh, prop.thresh, regression, low.threshold, up.threshold)
      ps_trimmed = microDecon_2_phyloseq(ps_obj = ps, env = theseus::pstoveg_sample(ps), decontaminated = decontaminated_ext, taxo_ranks = colnames(Taxo)[(colnames(Taxo)!="OTU_ID")])
    }else{
      decontaminated_ext <- microDecon::decon(data = ASV_table_ps, numb.blanks=sum(df_temp[df_temp$amplicon_type!="sample",]$n), numb.ind = df_temp[df_temp$amplicon_type=="sample",]$n, taxa = FALSE,runs, thresh, prop.thresh, regression, low.threshold, up.threshold)
      ps_trimmed = microDecon_2_phyloseq(ps_obj = ps, env = theseus::pstoveg_sample(ps), decontaminated = decontaminated_ext)
    }

  }else if(method == "max_v"){
    ps_blank = ps %>% subset_samples(amplicon_type!="sample")
    Extraction_neg_max_vec <- apply(ps_blank %>% theseus::pstoveg_otu() %>% t() %>% as.data.frame(), 1, max) %>% as.vector()
    names(Extraction_neg_max_vec) = taxa_names(ps_blank)
    Extractiondf = ps %>% theseus::pstoveg_otu() %>% as.data.frame() %>% dplyr::select(ASVs_in_Blanks)
    Extractiondf = sweep(Extractiondf,MARGIN=2,Extraction_neg_max_vec,FUN="-")
    Extractiondf <- replace(Extractiondf, Extractiondf < 0, 0)
    new_df = ps %>% theseus::pstoveg_otu() %>% as.data.frame() %>% dplyr::select(!ASVs_in_Blanks)
    new_df = cbind(new_df, Extractiondf)
    ps_trimmed = ps
    ps_trimmed@otu_table = otu_table(new_df, taxa_are_rows=FALSE)
    ps_trimmed = ps_trimmed %>% subset_samples(amplicon_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
  } else{
    ps_trimmed = ps %>% subset_samples(amplicon_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)
    ps_trimmed = ps_trimmed %>% pop_taxa(ASVs_in_Blanks)
  }
  # Printing results
  ntaxa_before = ps %>% subset_samples(amplicon_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% ntaxa()
  ntaxa_after = ntaxa(ps_trimmed)

  cat("\n")
  cat(paste0("Contamination removal outcome using ",method),"\n")
  cat(paste0("Number of ASVs totally removed: ",ntaxa_before - ntaxa_after))
  cat("\n")
  cat(paste0("Percent of ASVs removed: ",round((1 - (ntaxa_after / ntaxa_before)) * 100,2), " %"))
  cat("\n")

  reads_before = ps %>% subset_samples(amplicon_type == "sample") %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% sample_sums() %>% sum()
  reads_after = ps_trimmed %>% phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% sample_sums() %>% sum()

  cat(paste0("Total number of reads  removed: ",reads_before - reads_after))
  cat("\nPercent of reads removed: ",paste0(round((1 - (reads_after / reads_before)) * 100,2)," %"))
  cat("\n")
  return(ps_trimmed)
}
