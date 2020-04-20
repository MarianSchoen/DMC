#' deconvolute given bulks with DeconRNASeq using single cell data
#'
#' @param exprs numeric matrix containing single cell profiles as columns
#' @param pheno phenotype data corresponding to the expression matrix.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature list of cell types not to be included in the
#'  signature matrix
#' @param max.genes maximum number of genes that will be included in the
#'  signature for each celltype
#' @param optimize boolean, should the signature matrix be optimized by
#'  condition number? If FALSE, max.genes genes will be used
#' @param split.data boolean, should the training data be split for signature
#'  matrix creation? If TRUE, 10% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal features
#' @return list with one entry: est.props, matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' @example run_deconrnaseq(training.exprs, training.pheno, bulk.exprs)
run_deconrnaseq <- function(exprs, pheno, bulks, exclude.from.signature = NULL, max.genes = 500, optimize = TRUE, split.data = TRUE,
                            model = NULL) {
  suppressMessages(library(DeconRNASeq, quietly = TRUE))
  # error checking
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  if (nrow(exprs) != nrow(bulks)) {
      features <- intersect(rownames(exprs), rownames(bulks))
      if (length(features) > 0) {
          exprs <- exprs[features, ]
          bulks <- bulks[features, ]
      }
  }
  if (!is.null(max.genes) && max.genes == 0) {
      max.genes <- NULL
  }
  rownames(pheno) <- colnames(exprs)
  # scale to counts and create signature
  exprs <- scale_to_count(exprs)
  
  valid.model <- T
  if(!is.null(model)){
    if(all(c("ref.profiles", "g") %in% names(model))){
      full.mat <- model$ref.profiles
      g <- model$g
      if(all(names(g) %in% rownames(full.mat)) && length(g) == nrow(full.mat)){
        g <- g[rownames(full.mat)]
        if(any(duplicated(rownames(full.mat)))){
          g <- g[-which(duplicated(rownames(full.mat)))]
          full.mat <- full.mat[-which(duplicated(rownames(full.mat))),]
        }
        ref.profiles <- apply(full.mat, 2, function(x){x*g})
        if(any(rowSums(ref.profiles) == 0)){
          ref.profiles <- ref.profiles[-which(rowSums(ref.profiles) == 0),]
        }
        ref.profiles <- as.data.frame(ref.profiles)
      }else{
        warning("reference profiles and g vector do not contain the same genes")
        valid.model <- F
      }
    }else{
      warning("passed model parameter does not contain entries 'ref.profiles' and 'g'")
      valid.model <- F
    }
  }else{
    valid.model <- F
  }
  if(!valid.model){
    model <- create_sig_matrix(exprs,
        pheno,
        exclude.from.signature,
        max.genes = max.genes,
        optimize = optimize,
        split.data = split.data
      )
    
    ref.profiles <- as.data.frame(model$sig.matrix)
    full.mat <- model$full.matrix
    g <- model$g
  }
  # DeconRNASeq requires data frame as input
  df.mix <- as.data.frame(bulks)
  rownames(df.mix) <- rownames(bulks)

  # there is no option to switch the output of this function off...
  sink("/dev/null")
  result <- try(DeconRNASeq::DeconRNASeq(df.mix, ref.profiles))
  sink()
  if (!class(result) == "try-error") {
    result <- t(result$out.all[1:ncol(bulks), , drop = FALSE])
    colnames(result) <- colnames(bulks)
    if(!all(colnames(ref.profiles) %in% rownames(result))){
      result <- complete_estimates(result, colnames(ref.profiles))
    }
    return(list(est.props = result, sig.matrix = as.matrix(ref.profiles), ref.profiles = full.mat, g = g))
  } else {
    return(list(est.props = NULL, sig.matrix = as.matrix(ref.profiles), ref.profiles = full.mat, g = g))
  }
}
