#' deconvolute given bulks with BSEQ-sc using single cell data
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
#' @return list with one entry: est.props, matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' @example 
#' TODO: where to the 'training.exprs', 'training.pheno' and 
#' 'bulks.exprs' come from? 
#' 
#' run_bseqsc(
#'    training.exprs
#'    , training.pheno
#'    , bulks.exprs
#'  )
run_bseqsc <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = NULL,
  optimize = TRUE,
  split.data = TRUE
  ) {
  if(!exists("verbose")) verbose <- FALSE
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

  # do not normalize to counts for BSEQ-sc
  # exprs <- scale_to_count(exprs)

  if (!is.null(exclude.from.signature)) {
    to.exclude <- which(pheno[, "cell_type"] %in% exclude.from.signature)
    if  (length(to.exclude) > 0)  {
      exprs <- exprs[, -to.exclude]
      pheno <- pheno[-to.exclude, ]
    }
  }

  # use all supplied genes
  n.genes <- nrow(exprs)
  if (verbose) print("creating marker gene lists...")
  if (split.data) {
    learning.cells <- c()
    for (t in unique(pheno[, "cell_type"])) {
      learning.cells <- c(
        learning.cells,
        sample(
          which(pheno[, "cell_type"] == t),
          floor(0.9 * length(which(pheno[, "cell_type"] == t))),
          replace = FALSE
        )
      )
    }

    deg.per.type <- try(marker_genes(exprs[, learning.cells, drop = F], pheno[learning.cells, , drop = F], NULL))
    exprs <- exprs[, -learning.cells, drop = F]
    pheno <- pheno[-learning.cells, , drop = F]
  } else {
    deg.per.type <- try(marker_genes(exprs, pheno, NULL))
  }
  if (class(deg.per.type) == "try-error") {
    warning("Error while trying to find marker genes")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  if (length(deg.per.type) == 0) {
    warning("No genes passed the BSEQ-sc criteria")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # create single cell and bulk expression set
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

  B <- try({
    bseqsc_basis(sc.exprs, deg.per.type,
      clusters = "cell_type", samples = "sample.name",
      ct.scale = TRUE
    )
  })
  if (class(B) == "try-error") {
    warning("BSEQ-sc signature matrix creation failed")
    # save information in case of error
    save(sc.exprs, deg.per.type, file = paste("../bseqsc-error_", Sys.time(), sep = ""))
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  fit <- try(bseqsc_proportions(bulks, B, log = F, verbose = verbose))
  if (class(fit) == "try-error") {
    warning("BSEQ-sc estimation failed")
    return(list(est.props = NULL, sig.matrix = NULL))
  }
  if(!all(names(deg.per.type) %in% rownames(fit$coefficients))){
    est.props <- complete_estimates(fit$coefficients, names(deg.per.type))
  }else{
    est.props <- fit$coefficients
  }
  return(list(est.props = est.props, sig.matrix = B))
}
