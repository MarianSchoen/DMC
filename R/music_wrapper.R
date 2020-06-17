#' deconvolute given bulks with MuSiC using single cell data
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
#' @param verbose boolean
#' @param cell.type.column string, which column of 'training.pheno'/'test.pheno'
#' holds the cell type information?
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @return list with four entries: 
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm (features x cell types)
#' @example run_music(training.exprs, training.pheno, bulks[,1:5])
#'
run_music <- function(exprs,
                      pheno,
                      bulks,
                      exclude.from.signature = NULL,
                      max.genes = NULL,
                      optimize = TRUE,
                      split.data = TRUE,
                      cell.type.column = "cell_type",
                      patient.column = NULL
                      ) {
  suppressMessages(library(xbioc, quietly = TRUE))
	# parameters checks
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  features <- intersect(rownames(exprs), rownames(bulks))
  if (length(features) > 0) {
      exprs <- exprs[features, ]
      bulks <- bulks[features, ]
  }
  if(!is.null(max.genes)){
    if (max.genes == 0) {
        max.genes <- NULL
    }
  }
  if(is.null(patient.column) | !patient.column %in% colnames(pheno)){
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # MuSiC uses all supplied genes
  n.genes <- nrow(exprs)

  # ExpressionSet creation may fail without this...
  rownames(exprs) <- make.names(rownames(exprs))
  rownames(bulks) <- make.names(rownames(bulks))

  # create ExpressionSets from exprs, pheno and bulks
  sc.exprs <- Biobase::ExpressionSet(
    assayData = exprs,
    phenoData = Biobase::AnnotatedDataFrame(pheno)
  )
  bulk.pheno <- data.frame(sample = colnames(bulks))

  # equality check of sample names may fail due to different attributes
  rownames(bulk.pheno) <- colnames(bulks)
  colnames(bulks) <- rownames(bulk.pheno)
  bulks <- Biobase::ExpressionSet(
    assayData = bulks,
    phenoData = Biobase::AnnotatedDataFrame(bulk.pheno)
  )

  # exclude samples of types contained in exclude.from.signature
  if (!is.null(exclude.from.signature)) {
    if (length(which(unique(pheno[, cell.type.column]) %in% exclude.from.signature)) > 0) {
      include.in.x <- unique(pheno[, cell.type.column])[-which(unique(pheno[, cell.type.column]) %in% exclude.from.signature)]
    }else{
      include.in.x <- unique(pheno[, cell.type.column])
    }
  } else {
    include.in.x <- unique(pheno[, cell.type.column])
  }

  # deconvolution
  est.prop.music <- try(MuSiC::music_prop(
    bulk.eset = bulks, sc.eset = sc.exprs, clusters = cell.type.column,
    samples = patient.column, select.ct = include.in.x, verbose = FALSE
  ))
  if(class(est.prop.music) == "try-error"){
	  return(list(est.props = NULL, sig.matrix = NULL))
  }

  est.props <- as.matrix(t(est.prop.music$Est.prop.weighted))

  # complete the estimation matrix in case of cell type dropouts
  if(!all(include.in.x %in% rownames(est.props))){
    est.props <- complete_estimates(est.props, include.in.x)
  }
  
  return(list(est.props = est.props, sig.matrix = NULL))
}
