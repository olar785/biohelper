#' @title Pairwise dissimilarity per group.
#'
#' @param physeq A phyloseq-class object
#' @param group Grouping variable name (contained in \code{\link[phyloseq]{sample_data}})
#' @param between_groups Logical, estimate between-group dissimilarity
#' @param method Distance/dissimilarity method (as character string; see \code{\link[phyloseq]{distanceMethodList}})
#' @param method_title Logical, add method name to the plot title
#' @param notch Logical, draw the notches at each side of the boxes
#' @param justDF  Logical, instead of returning a ggplot2-object, return the data.frame that was used to build the plot
#' @param ... Additional arguments may be passed to \code{\link[phyloseq]{distance}} function from phyloseq package

#' @details This script is an adaptation of the 'phyloseq_group_dissimilarity' function of the metagMisc R package to include jaccard nestedness (jne) and turnover (jtu) as dissimilarity methods. If you use it, please cite the original package (citation('metagMisc')).
#'          If the notches of two boxplots do not overlap this indicates that the two medians differ.
#'
#' @return ggplot2-object or data.frame (if justDF = TRUE).
#' @export
#'
#' @examples
#' ## Load and subset data
#' data(enterotype)
#' ent <- subset_samples(enterotype, Enterotype %in% c("1", "2"))
#'
#' ## Dissimilarity
#' pw_group_dissimilarity(ent, group = "Enterotype", method = "jaccard")


pw_group_dissimilarity = function (physeq, group = NULL, between_groups = TRUE, method = "bray",
                                        method_title = FALSE, notch = TRUE, justDF = FALSE, ...)
{
  split_phyloseq_by_variable <- function(physeq, variable, drop_zeroes = TRUE) {
    mtd <- as.data.frame(phyloseq::sample_data(physeq))
    sample_ids <- phyloseq::sample_names(physeq)
    groups <- unique(mtd[[variable]])
    stats::setNames(lapply(groups, function(value) {
      keep <- sample_ids[mtd[[variable]] == value]
      out <- phyloseq::prune_samples(keep, physeq)
      if (drop_zeroes) {
        out <- phyloseq::prune_taxa(phyloseq::taxa_sums(out) > 0, out)
      }
      out
    }), groups)
  }

  dist_to_list <- function(x, tri = FALSE) {
    mat <- as.matrix(x)
    idx <- if (tri) which(lower.tri(mat), arr.ind = TRUE) else which(!is.na(mat), arr.ind = TRUE)
    data.frame(
      row = rownames(mat)[idx[, 1]],
      col = colnames(mat)[idx[, 2]],
      value = mat[idx],
      stringsAsFactors = FALSE
    )
  }

  if (is.null(phyloseq::sample_data(physeq, errorIfNULL = TRUE))) {
    stop("Error: Sample data is missing in the phyloseq-object.\n")
  }
  if (is.null(group)) {
    stop("Error: groupping variable should be specified.\n")
  }
  mtd <- as.data.frame(phyloseq::sample_data(physeq))
  mtd$SampleNames <- phyloseq::sample_names(physeq)
  if (!group %in% colnames(mtd)) {
    stop("Error: Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  tabb <- table(mtd[, group])
  if (any(is.na(mtd[, group]))) {
    stop("Error: there are NA values in the grouping variable.\n")
  }
  if (length(tabb) == 1) {
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  if (length(tabb) > 1) {
    physeq_split <- split_phyloseq_by_variable(physeq, variable = group, drop_zeroes = TRUE)
    ################ HERE
    if(method == "jtu"){
      dd <- plyr::llply(.data = physeq_split, .fun = function(z, ...) {
        betapart::beta.pair(z %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jtu"]]
      }, ...)
    }else if(method == "jne"){
      dd <- plyr::llply(.data = physeq_split, .fun = function(z, ...) {
        betapart::beta.pair(z %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jne"]]
      }, ...)
    }else{
      dd <- plyr::llply(.data = physeq_split, .fun = function(z, ...) {
        phyloseq::distance(z, method = method, type = "samples", ...)
      }, ...)
    }
    ################
    ddm <- plyr::ldply(.data = dd, .fun = function(z) {
      data.frame(Dist = as.vector(z), Comparison = "within-group")
    }, .id = "Group")

    if (between_groups == TRUE) {
      ################ HERE
      if(method == "jtu"){
        ddf <- betapart::beta.pair(physeq %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jtu"]]
      }else if(method == "jne"){
        ddf <- betapart::beta.pair(physeq %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jne"]]
      }else{
        ddf <- phyloseq::distance(physeq, method = method, type = "samples", ...)
      }
      ################

      ddl <- dist_to_list(ddf, tri = FALSE)
      ddl <- ddl[-with(ddl, which(col == row)), ]
      ddl$ColGroup <- mtd[match(x = ddl$col, table = mtd$SampleNames),
                          group]
      ddl$RowGroup <- mtd[match(x = ddl$row, table = mtd$SampleNames),
                          group]
      ddl <- ddl[which(!ddl$ColGroup == ddl$RowGroup),
      ]
      grr <- t(apply(ddl[, c("ColGroup", "RowGroup")],
                     1, sort))
      ddl$Group <- plyr::aaply(.data = grr, .margins = 1, .fun = paste,
                         collapse = "-")
      ddm <- rbind(ddm, data.frame(Group = ddl$Group,
                                   Dist = ddl$value, Comparison = "between-groups"))
    }
    pp <- ggplot2::ggplot(data = ddm, ggplot2::aes(x = Group, y = Dist, fill = Group)) +
      ggplot2::geom_boxplot(size = 0.8, notch = notch)
  }
  if (length(tabb) == 1) {

    ################ HERE
    if(method == "jtu"){
      dd <- betapart::beta.pair(physeq %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jtu"]]
    }else if(method == "jne"){
      dd <- betapart::beta.pair(physeq %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jne"]]
    }else{
      dd <- phyloseq::distance(physeq, method = method, type = "samples", ...)
    }
    ################

    ddm <- data.frame(Dist = as.vector(dd), Group = names(tabb)[1])
    pp <- ggplot2::ggplot(ddm, ggplot2::aes(x = Group, y = Dist)) + ggplot2::geom_boxplot(size = 0.8,
                                                               notch = notch)
  }
  pp <- pp + ggplot2::xlab(group) + ggplot2::ylab(paste("Pairwise dissimilarity (", method, ")", sep = ""))
  if (method_title == TRUE) {
    pp <- pp + ggplot2::ggtitle(method)
  }
  if (justDF == FALSE) {
    return(pp)
  }
  if (justDF == TRUE) {
    return(ddm)
  }
}
