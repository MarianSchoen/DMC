#' deconvolute given bulk with DTD using single-cell data
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
#' @param max.genes numeric, maximum number of genes that will be included in
#' the signature for each celltype, default 1000
#' @param verbose boolean, default FALSE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @param scale.cpm boolean, scale single-cell profiles to CPM? default FALSE
#' @param model pre-trained model for DTD deconvolution
#' as returned by this wrapper, default NULL
#' @param model_exclude character vector, cell type(s) to exclude
#' from the supplied pre-trained model, default NULL
#' @return list with four entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained \cr
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' '(features x cell types) \cr
#' 3) model - the trained DTD model\cr
#' @example run_dtd(training.exprs, training.pheno, bulks)
#' @export

run_dtd <- function(
  exprs,
  pheno,
  bulks,
  cell.type.column = "cell_type",
  exclude.from.signature = NULL,
  max.genes = 1000,
  verbose = FALSE,
  patient.column = NULL,
  scale.cpm = FALSE,
  model = NULL,
  model_exclude = NULL
) {
  suppressMessages(suppressWarnings(library(Matrix, quietly = TRUE)))
  # error checking
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  features <- intersect(rownames(exprs), rownames(bulks))
  if (length(features) > 0) {
      exprs <- exprs[features, ]
      bulks <- bulks[features, ]
  } else {
      stop("no common features in bulks and expression data.")
  }
  if (!is.null(max.genes) && max.genes == 0) {
      max.genes <- NULL
  }
  if (scale.cpm) {
    # prepare phenotype data and cell types to use
    exprs <- scale_to_count(exprs)
  }
  if (!is.matrix(exprs)) {
	  if (! class(exprs) == "dgCMatrix") {
		  stop("exprs must be a matrix or sparse matrix (dgCMatrix)")
	  }
  }
  if(!is.matrix(bulks)){
	  if(! class(bulks) == "dgCMatrix"){
		    stop("bulks must be a matrix or sparse matrix (dgCMatrix)")
	  }
  }

  # prepare phenotype data and cell types to use
  cell.types <- as.character(pheno[, cell.type.column])
  names(cell.types) <- colnames(exprs)

  cts <- unique(cell.types)
  include.in.x <- cts
  if (! is.null(exclude.from.signature)) {
    if (length(which(cts %in% exclude.from.signature)) > 0) {
      include.in.x <- cts[-which(cts %in% exclude.from.signature)]
    }
  }
  rm(cts)

  if(is.null(model)){
    if (is.null(max.genes)) max.genes <- 1000
    # create reference profiles
    sample.X <- DTD::sample_random_X(
      included.in.X = include.in.x,
      pheno = cell.types,
      expr.data = Matrix::as.matrix(exprs),
      percentage.of.all.cells = 0.3,
      normalize.to.count = TRUE
    )
    sig.matrix <- sample.X$X.matrix
    rm(sample.X)

    # do no remove used samples right now; try keeping reference samples
    # an option might be to remove these samples from the training data:
    # samples.to.remove <- sample.X$samples.to.remove

    # choose either max.genes genes per cell type or all available genes
    # but set maximum to 4000 due to runtime
    n.genes <- min(4000, nrow(exprs), length(unique(cell.types)) * max.genes)
    # select the top n.genes most variable genes for deconvolution
    top.features <- rownames(exprs)[
      order(apply(exprs, 1, var), decreasing = TRUE)[1:n.genes]
    ]
    exprs <- exprs[which(rownames(exprs) %in% top.features), ]
    sig.matrix <- sig.matrix[top.features, ]

    # create artificial mixtures to train DTD model on
    n.per.mixture <- max(floor(0.1 * ncol(exprs)),3)
    n.samples <- max(n.genes, 50)
    training.bulks <- DTD::mix_samples(
      expr.data = Matrix::as.matrix(exprs),
      pheno = cell.types,
      included.in.X = include.in.x,
      n.samples = n.samples,
      n.per.mixture = n.per.mixture,
      verbose = FALSE
    )

    # set the starting parameters to 1
    start.tweak <- rep(1, n.genes)
    names(start.tweak) <- top.features

    # train the DTD model
    dtd.model <- suppressMessages(try(
      DTD::train_deconvolution_model(
        tweak = start.tweak,
        X.matrix = sig.matrix,
        train.data.list = training.bulks,
        estimate.c.type = "direct",
        verbose = FALSE,
        NORM.FUN = "identity",
        #learning.rate = 1,
        cv.verbose = FALSE
      ),
      silent = TRUE
    ))
  }else{
    dtd.model <- model
    if (!is.null(model_exclude)) {
      cts <- colnames(dtd.model$reference.X)
      if (all(model_exclude %in% cts)) {
        cts <- cts[-which(cts %in% model_exclude)]
        dtd.model$reference.X <- dtd.model$reference.X[, cts, drop = FALSE]
      }else{
        stop("Not all cell types in 'model_exclude' are present in the model")
      }
    }
    sig.matrix <- dtd.model$reference.X
  }

  if (!class(dtd.model) == "try-error") {
    # use the model to estimate the composition of the supplied bulks
    est.props <- DTD::estimate_c(
      new.data = Matrix::as.matrix(bulks[rownames(sig.matrix), , drop = F]),
      DTD.model = dtd.model,
      estimate.c.type = "direct"
    )

    # if any cell types dropped out during estimation complete the matrix
    if(!all(include.in.x %in% rownames(est.props))){
      est.props <- complete_estimates(est.props, include.in.x)
    }
  } else {
    # if the model building failed,
    est.props <- NULL
  }

  # return estimated proportions and the effective signature matrix used
  return(list(est.props = est.props,
          sig.matrix = sig.matrix,
          model = dtd.model))
}
