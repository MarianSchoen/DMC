# written by Tim Mitus

#' deconvolute given bulks with DeconRNASeq using single cell data
#'
#' @@param exprs matrix containing single cell profiles as columns
#' @param pheno phenotype data corresponding to the expression matrix.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature list of cell types not to be included in the signature matrix
#' @param max.genes maximum number of genes that will be included in the signature for each celltype
#' @param optimize logical, should the signature matrix be optimized by condition number? If FALSE, max.genes genes will be used
#' @param split.data logical, should the training data be split for signature matrix creation? If TRUE, 10% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal features
#' @return list with one entry: est.props, matrix containing for each bulk the estimated fractions of the cell types contained
#' @example run_deconrnaseq(training.exprs, training.pheno, bulk.exprs)


suppressMessages(library(DeconRNASeq, quietly = TRUE))

run_deconrnaseq <- function(exprs, pheno, bulks, exclude.from.signature = NULL, max.genes = 500, optimize = TRUE, split.data = TRUE) {
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

  # scale to counts and create signature
  exprs <- scale_to_count(exprs)

  ref.profiles <- as.data.frame(create_sig_matrix(exprs,
      pheno,
      exclude.from.signature,
      max.genes = max.genes,
      optimize = optimize,
      split.data = split.data
    ))
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
    return(list(est.props = result, sig.matrix = as.matrix(ref.profiles)))
  } else {
    # store information of the error
    save(df.mix,
      ref.profiles,
      result,
      file = paste("../error_deconrnaseq_", Sys.time(), ".rda", sep = "")
    )
    return(list(est.props = NULL, sig.matrix = as.matrix(ref.profiles)))
  }
}
