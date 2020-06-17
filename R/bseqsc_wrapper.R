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
#' matrix creation? If TRUE, 10\% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal
#' features
#' @param verbose boolean
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information?
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @return list with four entries: 
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm (features x cell types)

run_bseqsc <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = NULL,
  optimize = TRUE,
  split.data = TRUE,
  verbose = FALSE,
  cell.type.column = "cell_type",
  patient.column = NULL
  ) {
	suppressMessages(library(Biobase, quietly = T))
	suppressMessages(library(xbioc, quietly = T))
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
  if(is.null(patient.column) | ! patient.column %in% colnames(pheno)){
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # ExpressionSet creation may fail in some cases without this
  rownames(exprs) <- make.names(rownames(exprs))
  rownames(bulks) <- make.names(rownames(bulks))
  
  bulk.pheno <- data.frame(colnames(bulks))
  rownames(bulk.pheno) <- colnames(bulks)
  colnames(bulk.pheno) <- "sample"
  colnames(bulks) <- rownames(bulk.pheno)   # ExpressionSet creation may fail due to differences in attributes
  bulks <- Biobase::ExpressionSet(
    assayData = bulks,
    phenoData = Biobase::AnnotatedDataFrame(bulk.pheno)
  )

  # exclude cells of types contained in exclude.from.signature
  if (!is.null(exclude.from.signature)) {
    to.exclude <- which(pheno[, cell.type.column] %in% exclude.from.signature)
    if  (length(to.exclude) > 0)  {
      exprs <- exprs[, -to.exclude]
      pheno <- pheno[-to.exclude, ]
    }
  }

  # use all supplied genes
  n.genes <- nrow(exprs)
  if (verbose) print("creating marker gene lists...")
  # sample cells to use for gene selection in case of split training data and
  # find informative genes for each cell type
  if (split.data) {
    learning.cells <- c()
    for (t in unique(pheno[, cell.type.column])) {
      learning.cells <- c(
        learning.cells,
        sample(
          which(pheno[, cell.type.column] == t),
          floor(0.7 * length(which(pheno[, cell.type.column] == t))),
          replace = FALSE
        )
      )
    }

    deg.per.type <- try(marker_genes(exprs[, learning.cells, drop = F], pheno[learning.cells, , drop = F], NULL, cell.type.column = cell.type.column))
    exprs <- exprs[, -learning.cells, drop = F]
    pheno <- pheno[-learning.cells, , drop = F]
  } else {
    deg.per.type <- try(marker_genes(exprs, pheno, NULL, cell.type.column = cell.type.column))
  }

  # BSEQ-sc marker selection produces errors sometimes
  if (class(deg.per.type) == "try-error") {
    warning("Error while trying to find marker genes for BSEQ-sc")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  if (length(deg.per.type) == 0) {
    warning("No genes passed the BSEQ-sc criteria")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # create single cell and bulk expression set
  sc.exprs <- Biobase::ExpressionSet(
    assayData = exprs,
    phenoData = Biobase::AnnotatedDataFrame(pheno)
  )

  # create reference matrix based on single cell data and DEGs
  B <- try({
	  bseqsc::bseqsc_basis(sc.exprs, deg.per.type,
      clusters = cell.type.column, samples = patient.column,
      ct.scale = TRUE
    )
  })
  
  if (class(B) == "try-error") {
    warning("BSEQ-sc signature matrix creation failed")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # deconvolute bulks using reference matrix B
  fit <- try(bseqsc::bseqsc_proportions(bulks, B, log = F, verbose = verbose))

  if (class(fit) == "try-error") {
    warning("BSEQ-sc estimation failed")
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # complete the estimation matrix in case of dropout cell types
  if(!all(unique(pheno[[cell.type.column]]) %in% rownames(fit$coefficients))){
    est.props <- bseqsc::complete_estimates(fit$coefficients, unique(pheno[[cell.type.column]]))
  }else{
    est.props <- fit$coefficients
  }
  return(list(est.props = est.props, sig.matrix = B))
}
