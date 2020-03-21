# written by Tim Mirus

#' deconvolute given bulks with MuSiC using single cell data
#'
#' @param exprs matrix containing single cell profiles as columns
#' @param pheno phenotype data corresponding to the expression matrix.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature list of cell types not to be included in the signature matrix
#' @param max.genes maximum number of genes that will be included in the signature for each celltype
#' @param optimize logical, should the signature matrix be optimized by condition number? If FALSE, max.genes genes will be used
#' @param split.data logical, should the training data be split for signature matrix creation? If TRUE, 10% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal features
#' @return list with one entry: est.props, matrix containing for each bulk the estimated fractions of the cell types contained
#' @example run_music(training.exprs, training.pheno, bulks[,1:5])
#'
run_music <- function(exprs,
                      pheno,
                      bulks,
                      exclude.from.signature = NULL,
                      max.genes = NULL,
                      optimize = TRUE,
                      split.data = TRUE) {
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
  if(!is.null(max.genes)){
    if (max.genes == 0) {
        max.genes <- NULL
    }
  }

  # MuSiC uses all supplied genes
  n.genes <- nrow(exprs)

  # do not normalize profiles to counts
  # exprs <- scale_to_count(exprs)

  # create ExpressionSets from exprs, pheno and bulks
  sc.exprs <- ExpressionSet(
    assayData = exprs,
    phenoData = AnnotatedDataFrame(pheno)
  )
  bulk.pheno <- data.frame(colnames(bulks))
  rownames(bulk.pheno) <- colnames(bulks)
  colnames(bulk.pheno) <- "sample"

  bulks <- ExpressionSet(
    assayData = bulks,
    phenoData = AnnotatedDataFrame(bulk.pheno)
  )
  if (!is.null(exclude.from.signature)) {
    if (length(which(unique(pheno[, "cell_type"]) %in% exclude.from.signature)) > 0) {
      include.in.x <- unique(pheno[, "cell_type"])[-which(unique(pheno[, "cell_type"]) %in% exclude.from.signature)]
    }else{
      include.in.x <- unique(pheno[, "cell_type"])
    }
  } else {
    include.in.x <- unique(pheno[, "cell_type"])
  }

  # deconvolution
  est.prop.music <- try(music_prop(
    bulk.eset = bulks, sc.eset = sc.exprs, clusters = "cell_type",
    samples = "sample.name", select.ct = include.in.x, verbose = FALSE
  ))
  if(class(est.prop.music) == "try-error"){
	  return(list(est.props = NULL, sig.matrix = NULL))
  }

  est.props <- as.matrix(t(est.prop.music$Est.prop.weighted))
  return(list(est.props = est.props, sig.matrix = NULL))
}
