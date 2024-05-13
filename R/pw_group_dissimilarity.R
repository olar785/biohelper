#' @title Pairwise dissimilarity per group.
#'
#' @param physeq A phyloseq-class object
#' @param group Grouping variable name (contained in \code{\link{sample_data}})
#' @param between_groups Logical, estimate between-group dissimilarity
#' @param method Distance/dissimilarity method (as character string; see \code{\link{distanceMethodList}})
#' @param method_title Logical, add method name to the plot title
#' @param notch Logical, draw the notches at each side of the boxes
#' @param justDF  Logical, instead of returning a ggplot2-object, return the data.frame that was used to build the plot
#' @param ... Additional arguments may be passed to \code{\link{distance}} function from phyloseq package

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
#' ## Dissimilarity boxplots
#' phyloseq_group_dissimilarity(ent, group = "Enterotype")
#' phyloseq_group_dissimilarity(ent, group = "Enterotype", between_groups = F)
#' phyloseq_group_dissimilarity(ent, group = "Enterotype", method = "jaccard")


pw_group_dissimilarity = function (physeq, group = NULL, between_groups = TRUE, method = "bray",
                                        method_title = FALSE, notch = TRUE, justDF = FALSE, ...)
{
  require(plyr)
  if (justDF == FALSE) {
    require(ggplot2)
  }
  if (is.null(sample_data(physeq, errorIfNULL = T))) {
    stop("Error: Sample data is missing in the phyloseq-object.\n")
  }
  if (is.null(group)) {
    stop("Error: groupping variable should be specified.\n")
  }
  mtd <- as(object = sample_data(physeq), Class = "data.frame")
  mtd$SampleNames <- sample_names(physeq)
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
    physeq_split <- phyloseq_sep_variable(physeq, variable = group, drop_zeroes = T)
    ################ HERE
    if(method == "jtu"){
      dd <- llply(.data = physeq_split, .fun = function(z, ...) {
        betapart::beta.pair(z %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jtu"]]
      }, ...)
    }else if(method == "jne"){
      dd <- llply(.data = physeq_split, .fun = function(z, ...) {
        betapart::beta.pair(z %>% microbiome::transform("pa") %>% pstoveg_otu(), index.family = "jaccard",...)[["beta.jne"]]
      }, ...)
    }else{
      dd <- llply(.data = physeq_split, .fun = function(z, ...) {
        phyloseq::distance(z, method = method, type = "samples", ...)
      }, ...)
    }
    ################
    ddm <- ldply(.data = dd, .fun = function(z) {
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

      ddl <- dist2list(ddf, tri = FALSE)
      ddl <- ddl[-with(ddl, which(col == row)), ]
      ddl$ColGroup <- mtd[match(x = ddl$col, table = mtd$SampleNames),
                          group]
      ddl$RowGroup <- mtd[match(x = ddl$row, table = mtd$SampleNames),
                          group]
      ddl <- ddl[which(!ddl$ColGroup == ddl$RowGroup),
      ]
      grr <- t(apply(ddl[, c("ColGroup", "RowGroup")],
                     1, sort))
      ddl$Group <- aaply(.data = grr, .margins = 1, .fun = paste,
                         collapse = "-")
      ddm <- rbind(ddm, data.frame(Group = ddl$Group,
                                   Dist = ddl$value, Comparison = "between-groups"))
    }
    pp <- ggplot(data = ddm, aes(x = Group, y = Dist, fill = Group)) +
      geom_boxplot(size = 0.8, notch = notch)
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
    pp <- ggplot(ddm, aes(x = Group, y = Dist)) + geom_boxplot(size = 0.8,
                                                               notch = notch)
  }
  pp <- pp + xlab(group) + ylab(paste("Pairwise dissimilarity (", method, ")", sep = ""))
  if (method_title == TRUE) {
    pp <- pp + ggtitle(method)
  }
  if (justDF == FALSE) {
    return(pp)
  }
  if (justDF == TRUE) {
    return(ddm)
  }
}
