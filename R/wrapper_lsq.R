#' deconvolute given bulk with DTD using single-cell data
#' without loss function learning
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
#' @param max.genes numeric, maximum number of genes that will be included in
#' the signature for each celltype, default 500
#' @param verbose boolean, default FALSE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @param scale.cpm boolean, scale single-cell profiles to CPM? default FALSE
#' @param model pre-trained model for LSQ deconvolution
#' as returned by this wrapper, default NULL
#' @param model_exclude character vector, cell type(s) to exclude
#' from the supplied pre-trained model, default NULL
#' @return list with three entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained\cr
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' (features x cell types)\cr
#' 3) model - list containing reference.X (signature matrix)
#' and g (weight vector, all weights 1)\cr
#' @example run_dtd_baseline(training.exprs, training.pheno, bulks)
#' @export
run_least_squares <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = 500,
  cell.type.column = "cell_type",
  patient.column = NULL,
  scale.cpm = FALSE,
  model = NULL,
  model_exclude = NULL
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
  if (!is.null(max.genes)) {
    if (max.genes == 0) {
        max.genes <- NULL
    }
  }

  cell.types <- as.character(pheno[, cell.type.column])
  names(cell.types) <- colnames(exprs)

  cts <- unique(cell.types)
  include.in.x <- cts
  # exclude samples of types contained in exclude.from.signature
  if (!is.null(exclude.from.signature)) {
    if (length(which(cts %in% exclude.from.signature)) > 0) {
      include.in.x <- cts[-which(cts %in% exclude.from.signature)]
    }
  }

  if (is.null(model)) {
    if(scale.cpm){
      # prepare phenotype data and cell types to use
      exprs <- scale_to_count(exprs)
    }
    # create reference profiles
    sample.X <- DTD::sample_random_X(
      included.in.X = include.in.x,
      pheno = cell.types,
      expr.data = Matrix::as.matrix(exprs),
      percentage.of.all.cells = 0.9999,
      normalize.to.count = TRUE
    )
    sig.matrix <- sample.X$X.matrix
    rm(sample.X)
    # do not remove used samples from expression matrix and pheno data

    # use a maximum of 4000 genes
    n.genes <- min(4000, nrow(exprs))
    # use the top n.genes most variable genes
    top.features <- rownames(exprs)[
      order(apply(exprs, 1, var), decreasing = TRUE)[1:n.genes]
    ]
    exprs <- exprs[which(rownames(exprs) %in% top.features), ]
    sig.matrix <- sig.matrix[top.features, ]

    # set the model without learning
    start.tweak <- rep(1, n.genes)
    names(start.tweak) <- top.features

    model <- list(reference.X = sig.matrix, g = start.tweak)
  }else{
    sig.matrix <- model$reference.X
    start.tweak <- model$g
    if (!is.null(model_exclude)) {
      cts <- colnames(sig.matrix)
      if (all(model_exclude %in% cts)) {
        cts <- cts[-which(cts %in% model_exclude)]
        sig.matrix <- sig.matrix[, cts, drop = FALSE]
      }else{
        stop("Not all cell types in 'model_exclude' are present in the model")
      }
    }
  }

  # use the untrained model to estimate the composition of the supplied bulks
  est.props <- DTD::estimate_c(
    X.matrix = sig.matrix,
    new.data = Matrix::as.matrix(bulks[rownames(sig.matrix), , drop = FALSE]),
    DTD.model = start.tweak,
    estimate.c.type = "direct"
  )

  # if any cell types dropped out during estimation complete the matrix
  if (!all(include.in.x %in% rownames(est.props))) {
    est.props <- complete_estimates(est.props, include.in.x)
  }

  return(list(
    est.props = est.props,
    sig.matrix = sig.matrix,
    model = model
  ))
}
