#' deconvolute given bulk with DTD using single-cell data without loss function learning
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
#' @return list with four entries: 
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm (features x cell types)
#' @example run_dtd_baseline(training.exprs, training.pheno, bulks)
run_least_squares <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = 500,
  optimize = TRUE,
  split.data = TRUE, 
  cell.type.column = "cell_type"
) {
  # error checking
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
  # prepare phenotype data and cell types to use
  exprs <- scale_to_count(exprs)
  
  cell.types <- as.character(pheno[, cell.type.column])
  names(cell.types) <- colnames(exprs)

  # exclude samples of types contained in exclude.from.signature
  if (!is.null(exclude.from.signature)) {
    if (length(which(cell.types %in% exclude.from.signature)) > 0) {
      include.in.x <- unique(cell.types[-which(cell.types %in% exclude.from.signature)])
    }else{
      include.in.x <- unique(cell.types)
    }
  } else {
    include.in.x <- unique(cell.types)
  }
  
  # create reference profiles
  sample.X <- DTD::sample_random_X(
    included.in.X = include.in.x,
    pheno = cell.types,
    expr.data = exprs,
    percentage.of.all.cells = 0.3,
    normalize.to.count = TRUE
  )
  sig.matrix <- sample.X$X.matrix
  full.mat <- sig.matrix

  # remove used samples from expression matrix and pheno data
  #samples.to.remove <- sample.X$samples.to.remove
  #exprs <- exprs[, -which(colnames(exprs) %in% samples.to.remove)]
  #cell.types <- cell.types[-which(names(cell.types) %in% samples.to.remove)]

  # use a maximum of 4000 genes
  n.genes <- min(4000, nrow(exprs))
  # use the top n.genes most variable genes
  top.features <- rownames(exprs)[order(apply(exprs, 1, var), decreasing = TRUE)[1:n.genes]]
  exprs <- exprs[top.features, ]
  sig.matrix <- sig.matrix[top.features, ]

  # set the model without learning
  start.tweak <- rep(1, n.genes)
  names(start.tweak) <- top.features
    
  # use the untrained model to estimate the composition of the supplied bulks
  est.props <- DTD::estimate_c(
    X.matrix = sig.matrix,
    new.data = bulks[rownames(sig.matrix), , drop = FALSE],
    DTD.model = start.tweak,
    estimate.c.type = "direct"
  )

  # if any cell types dropped out during estimation complete the matrix
  if(!all(include.in.x %in% rownames(est.props))){
    est.props <- complete_estimates(est.props, include.in.x)
  }
  
  return(list(
    est.props = est.props,
    sig.matrix = sig.matrix
  ))
}
