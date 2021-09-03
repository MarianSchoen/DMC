#' deconvolute given bulks with MuSiC using single cell data
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
#' @param max.genes unused, MuSiC uses all genes
#' @param verbose boolean, default FALSE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default "patient"
#' @param scale.cpm boolean, scale single-cell profiles to CPM? default FALSE
#' @param model pre-trained model for MuSiC deconvolution
#' as returned by this wrapper, default NULL
#' @param model_exclude character vector, cell type(s) to exclude
#' from the supplied pre-trained model, default NULL
#' @return list with three entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained\cr
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' (features x cell types)\cr
#' 3) model - list containing learned parameters for music.basic function call
#' @example run_music(training.exprs, training.pheno, bulks[,1:5])
#' @export
run_music <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = NULL,
  cell.type.column = "cell_type",
  patient.column = "patient",
  scale.cpm = FALSE,
  model = NULL,
  model_exclude = NULL
  ) {
  suppressWarnings(suppressMessages(library(xbioc, quietly = TRUE)))
	# parameters checks
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  features <- intersect(rownames(exprs), rownames(bulks))
  if (length(features) > 0) {
      exprs <- exprs[features, ]
      bulks <- bulks[features, ]
  }
  if (is.null(patient.column)) {
	  cat("Patient column variable not present\n")
    return(list(est.props = NULL, sig.matrix = NULL))
  }
  if (!patient.column %in% colnames(pheno)) {
	  cat("Patient column not present\n")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  if (scale.cpm) {
    # prepare phenotype data and cell types to use
    exprs <- scale_to_count(exprs)
  }
  # MuSiC uses all supplied genes
  exprs <- Matrix::as.matrix(exprs)
  bulks <- Matrix::as.matrix(bulks)
  colnames(bulks) <- make.names(colnames(bulks))

  # ExpressionSet creation may fail without this...
  colnames(exprs) <- make.names(colnames(exprs))
  rownames(exprs) <- make.names(rownames(exprs))
  rownames(bulks) <- make.names(rownames(bulks))
  rownames(pheno) <- colnames(exprs)
  Y <- bulks

  cts <- unique(pheno[, cell.type.column])
  include.include.in.x <- cts
  # exclude samples of types contained in exclude.from.signature
  if (!is.null(exclude.from.signature)) {
    if (length(which(cts %in% exclude.from.signature)) > 0) {
      include.in.x <- cts[-which(cts %in% exclude.from.signature)]
    }
  }

  if (is.null(model)) {
    # create ExpressionSets from exprs, pheno and bulks
    sc.exprs <- Biobase::ExpressionSet(
      assayData = exprs,
      phenoData = Biobase::AnnotatedDataFrame(pheno)
    )
    bulk.pheno <- data.frame(sample = colnames(bulks))
    rm(exprs)

    # equality check of sample names may fail due to different attributes
    rownames(bulk.pheno) <- colnames(bulks)
    colnames(bulks) <- rownames(bulk.pheno)


    if(any(table(pheno[, cell.type.column])[include.in.x] < 2)){
      # remove cell types with only one sample
      include.in.x <- include.in.x[
        -which(table(pheno[, cell.type.column])[include.in.x] < 2)
      ]
    }
    basis <- MuSiC::music_basis(
      x = sc.exprs, clusters = cell.type.column,
      samples = patient.column, select.ct = include.in.x, verbose = FALSE
    )

    X <- basis$Disgn.mtx
    mS <- basis$M.S
    sigma <- basis$Sigma
    iter.max <- 1000
    nu <- 1e-04
    eps <- 0.01

    genes <- rownames(sigma)[
      which(!apply(sigma, 1, function(x){any(is.na(x))}))
    ]
    sigma <- sigma[genes,]
    X <- X[genes,]
    Y <- Y[genes,]

    # save model for later use
    model <- list(X = X, mS = mS, sigma = sigma, iter.max = iter.max, nu = nu, eps = eps)
  }else{
    # load from model variable
    X <- model$X
    mS <- model$mS
    sigma <- model$sigma
    iter.max <- model$iter.max
    nu <- model$nu
    eps <- model$eps
    Y <- Y[rownames(X),]

    if (!is.null(model_exclude)) {
      if (all(model_exclude %in% colnames(X))) {
        cts <- colnames(X)[-which(colnames(X) %in% model_exclude)]
        X <- X[, cts, drop = FALSE]
        mS <- mS[cts]
        sigma <- sigma[, cts, drop = FALSE]
      }else{
        stop("Not all cell types in 'model_exclude' are present in the model")
      }
    }
  }

  # deconvolution
  est.props <- apply(Y, 2, function(y){
    MuSiC::music.basic(
      Y = y,
      X = X,
      S = mS,
      Sigma = sigma,
      iter.max = iter.max,
      nu = nu,
      eps = eps
    )$p.weight
  })
  rownames(est.props) <- colnames(X)

  # complete the estimation matrix in case of cell type dropouts
  if(!all(include.in.x %in% rownames(est.props))){
    est.props <- complete_estimates(est.props, include.in.x)
  }

  return(list(
    est.props = est.props,
    sig.matrix = NULL,
    model = model)
  )
}
