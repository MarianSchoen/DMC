#' deconvolute given bulk with DTD using single-cell data
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'. 
#' Has to contain single cell labels in 'cell.type.column'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
#' @param max.genes numeric, maximum number of genes that will be included in 
#' the signature for each celltype
#' @param optimize boolean, should the signature matrix be optimized by
#' condition number? If FALSE, max.genes genes will be used
#' @param split.data boolean, should the training data be split for signature
#' matrix creation? If TRUE, 10\% of the data will be used to build
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
#' @example run_dtd(training.exprs, training.pheno, bulks)
#' @export

run_dtd <- function(
  exprs,
  pheno,
  bulks,
  cell.type.column = "cell_type", 
  exclude.from.signature = NULL,
  max.genes = NULL,
  optimize = TRUE,
  split.data = TRUE,
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
  if(scale.cpm){
    # prepare phenotype data and cell types to use
    exprs <- scale_to_count(exprs)
  }
  if(!is.matrix(exprs)){
	  if(class(exprs) == "dgCMatrix"){
		  #exprs <- as.matrix(exprs)
	  }else{
		  stop("exprs must be a matrix or sparse matrix (dgCMatrix)")
	  }
  }
  if(!is.matrix(bulks)){
	  if(class(bulks) == "dgCMatrix"){
		#bulks <- as.matrix(bulks)
	  }else{
		stop("bulks must be a matrix or sparse matrix (dgCMatrix)")
	  }
  }
  
  
    # prepare phenotype data and cell types to use
    cell.types <- as.character(pheno[, cell.type.column])
    names(cell.types) <- colnames(exprs)
    
    if (!is.null(exclude.from.signature)) {
      if (length(which(cell.types %in% exclude.from.signature)) > 0) {
        include.in.x <- unique(cell.types[-which(cell.types %in% exclude.from.signature)])
      }else{
        include.in.x <- unique(cell.types)
      }
    } else {
      include.in.x <- unique(cell.types)
    }
  if(is.null(model)){
    if (is.null(max.genes)) max.genes <- 4000
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
  
    # remove used samples from expression matrix and pheno data
    # do no remove samples right now; try keeping reference samples
    # samples.to.remove <- sample.X$samples.to.remove
  
    # # remove samples only if there is more than one of this cell type
    # if (any(table(cell.types) == 1)) {
    #   samples.to.retain <- c()
    #   for (t in include.in.x) {
    #     if (table(cell.types[samples.to.remove])[t] == 1) {
    #       samples.to.retain <- c(
    #         samples.to.retain,
    #         which(cell.types[samples.to.remove] == t)
    #       )
    #     }
    #   }
    #   samples.to.remove <- samples.to.remove[-samples.to.retain]
    # }
    # if(any(colnames(exprs) %in% samples.to.remove)){
    #   exprs <- exprs[, -which(colnames(exprs) %in% samples.to.remove)]
    # }
    # if(any(names(cell.types) %in% samples.to.remove)){
    #   cell.types <- cell.types[-which(names(cell.types) %in% samples.to.remove)]
    # }
  
    # choose either max.genes genes per cell type or all available genes
    # but set maximum to 4000 due to runtime
    n.genes <- min(4000, nrow(exprs), length(unique(cell.types)) * max.genes)
    # select the top n.genes most variable genes for deconvolution
    top.features <- rownames(exprs)[order(apply(exprs, 1, var), decreasing = TRUE)[1:n.genes]]
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
    dtd.model <- suppressMessages(try(DTD::train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = sig.matrix,
      train.data.list = training.bulks,
      estimate.c.type = "direct",
      verbose = FALSE,
      NORM.FUN = "identity",
      #learning.rate = 1,
      cv.verbose = FALSE
    ), silent = TRUE))
  }else{
    dtd.model <- model
    if(!is.null(model_exclude)){
      if(all(model_exclude %in% colnames(dtd.model$reference.X))){
        dtd.model$reference.X <- dtd.model$reference.X[,-which(colnames(dtd.model$reference.X) %in% model_exclude), drop = FALSE]
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
    #if(!all(include.in.x %in% rownames(est.props))){
    #  est.props <- complete_estimates(est.props, include.in.x)
    #}
  } else {
    # if the model building failed, 
    est.props <- NULL
  }

  # return estimated proportions and the effective signature matrix used
  return(list(est.props = est.props,
          sig.matrix = sig.matrix,
          model = dtd.model))
}
