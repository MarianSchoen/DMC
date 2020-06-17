#' deconvolute given bulks with DeconRNASeq using single cell data
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'. 
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
#' @param max.genes numeric, maximum number of genes that will be included in 
#' the signature for each celltype
#' @param optimize boolean, should the signature matrix be optimized by
#' condition number? If FALSE, max.genes genes will be used
#' @param split.data boolean, should the training data be split for signature
#' matrix creation? If TRUE, 10% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal
#' features
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? 
#' @return list with two entries: 
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm (features x cell types)
#' @example run_deconrnaseq(training.exprs, training.pheno, bulk.exprs)
run_deconrnaseq <- function(exprs, pheno, bulks, exclude.from.signature = NULL, max.genes = 500, optimize = TRUE, split.data = TRUE, cell.type.column = "cell_type") {
  # error checking
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  features <- intersect(rownames(exprs), rownames(bulks))
  if (length(features) > 0) {
      exprs <- exprs[features, ]
      bulks <- bulks[features, ]
  }
  if (!is.null(max.genes) && max.genes == 0) {
      max.genes <- NULL
  }

  # scale to counts
  exprs <- scale_to_count(exprs)
  
  # create signature matrix (DeconRNASeq needs data frames)
  ref.profiles <- create_sig_matrix(exprs,
      pheno,
      exclude.from.signature,
      max.genes = max.genes,
      optimize = optimize,
      split.data = split.data,
      cell.type.column = cell.type.column
  )

  # create bulk data frame
  df.mix <- as.data.frame(bulks)
  rownames(df.mix) <- rownames(bulks)

  # there is no option to switch the output of this function off
  # deconvolute
  invisible(capture.output(result <- try(DeconRNASeq::DeconRNASeq(df.mix, as.data.frame(ref.profiles)), silent = TRUE)))

  if (!class(result) == "try-error") {
    # select the interesting rows and rotate to be compatible with other algorithms' outputs
    result <- t(result$out.all[1:ncol(bulks), , drop = FALSE])
    colnames(result) <- colnames(bulks)

    # complete estimation matrix in case of droput cell types
    if(!all(colnames(ref.profiles) %in% rownames(result))){
      result <- complete_estimates(result, colnames(ref.profiles))
    }
    return(list(est.props = result, sig.matrix = ref.profiles))
  } else {
    return(list(est.props = NULL, sig.matrix = ref.profiles))
  }
}
