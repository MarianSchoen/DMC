#' deconvolute given bulk with DTD using single-cell data without loss function learning
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
#' @return list containing 3 elements:
#' 1) est.props - matrix containing the estimated proportions of all cell types in each bulk as returned by dtd (cell type x bulk)
#' 2) est.props.colscake - est.props with columns rescaled to sum up to 1
#' 3) est.props.rowscale - est.props with rows rescaled so that maximum in each row is 1
#' @example run_dtd_baseline(training.exprs, training.pheno, bulks)

suppressMessages(library(DTD))

run_least_squares <- function(exprs,
                             pheno,
                             bulks,
                             exclude.from.signature = NULL,
                             max.genes = 500,
                             optimize = TRUE,
                             split.data = TRUE,
                             model = NULL) {
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
  # prepare phenotype data and cell types to use
  exprs <- scale_to_count(exprs)
  
  cell.types <- as.character(pheno[, "cell_type"])
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
  
  valid.model <- T
  if(!is.null(model)){
    if(all(c("ref.profiles", "g") %in% names(model))){
      full.mat <- model$ref.profiles
      g <- model$g
      if(all(names(g) %in% rownames(full.mat)) && length(g) == nrow(full.mat)){
        g <- g[rownames(full.mat)]
        sig.matrix <- apply(full.mat, 2, function(x){x*g})
        if(any(rowSums(sig.matrix) == 0)){
          start.tweak <- g[-which(rowSums(sig.matrix) == 0)]
          sig.matrix <- sig.matrix[-which(rowSums(sig.matrix) == 0),]
        }
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
    # create reference profiles
    sample.X <- sample_random_X(
      included.in.X = include.in.x,
      pheno = cell.types,
      expr.data = exprs,
      percentage.of.all.cells = 0.1,
      normalize.to.count = TRUE
    )
    sig.matrix <- sample.X$X.matrix
    full.mat <- sig.matrix
  
    # remove used samples from expression matrix and pheno data
    samples.to.remove <- sample.X$samples.to.remove
    exprs <- exprs[, -which(colnames(exprs) %in% samples.to.remove)]
    cell.types <- cell.types[-which(names(cell.types) %in% samples.to.remove)]
  
    n.genes <- min(4000, nrow(exprs))
    top.features <- rownames(exprs)[order(apply(exprs, 1, var), decreasing = TRUE)[1:n.genes]]
    exprs <- exprs[top.features, ]
    sig.matrix <- sig.matrix[top.features, ]
  
    n.per.mixture <- floor(0.1 * ncol(exprs))
    n.samples <- max(n.genes, 50)
  
    # set the model without learning
    start.tweak <- rep(1, n.genes)
    names(start.tweak) <- top.features
  }
  # use the model to estimate the composition of the supplied bulks
  est.props <- estimate_c(
    X.matrix = sig.matrix,
    new.data = bulks[rownames(sig.matrix), , drop = FALSE],
    DTD.model = start.tweak,
    estimate.c.type = "direct"
  )
  # rescale by column and by row in order to determine how this changes results
  # est.props.colscale <- est.props
  # est.props.rowscale <- est.props
  # if (any(est.props < 0)) {
  #   est.props.rowscale[est.props.rowscale < 0] <- 0
  #   est.props.colscale[est.props.colscale < 0] <- 0
  # }
  # est.props.colscale <- apply(est.props.colscale, 2, function(x) {
  #   x / sum(x)
  # })
  # est.props.rowscale <- t(apply(est.props.rowscale, 1, function(x) {
  #   x / max(x)
  # }))
  # rownames(est.props.rowscale) <- rownames(est.props.colscale)
  # if(!all(include.in.x %in% rownames(est.props))){
  #   est.props <- complete_estimates(est.props, include.in.x)
  # }
  
  g <- rep(0, nrow(full.mat))
  names(g) <- rownames(full.mat)
  g[rownames(sig.matrix)] <- 1
  
  return(list(
    est.props = est.props,
    sig.matrix = sig.matrix, 
    ref.profiles = full.mat, 
    g = g
  ))
}
