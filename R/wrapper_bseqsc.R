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
#' the signature for each celltype, default 500
#' @param verbose boolean
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default "patient"
#' @param scale.cpm boolean, unused
#' @param model model for BSEQ-sc deconvolution as returned by this wrapper,
#' default NULL
#' @param model_exclude character vector, cell type(s) to exclude
#' from the supplied pre-trained model, default NULL
#' @return list with four entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained\cr
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' (features x cell types)\cr
#' 3) model - list containing signature matrix \cr
#' @export

run_bseqsc <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = NULL,
  verbose = FALSE,
  cell.type.column = "cell_type",
  patient.column = "patient",
  scale.cpm = FALSE,
  model = NULL,
  model_exclude = NULL
) {
  suppressMessages(library(Biobase, quietly = TRUE))
  suppressMessages(library(xbioc, quietly = TRUE))
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
  if (!is.null(max.genes)) {
    if (max.genes == 0) {
        max.genes <- NULL
    }
  }
  if (is.null(patient.column)) {
    warning("No patient column provided")
    return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
  }
  if (!patient.column %in% colnames(pheno)) {
    warning(paste0("column '", patient.column, "' not found in pheno data"))
    return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
  }

  # ExpressionSet creation may fail in some cases without this
  rownames(exprs) <- make.names(rownames(exprs))
  rownames(bulks) <- make.names(rownames(bulks))

  bulk.pheno <- data.frame(colnames(bulks))
  rownames(bulk.pheno) <- colnames(bulks)
  colnames(bulk.pheno) <- "sample"
  # ExpressionSet creation may fail due to differences in attributes
  colnames(bulks) <- rownames(bulk.pheno)
  # create bulk expression set
  bulks <- Biobase::ExpressionSet(
    assayData = Matrix::as.matrix(bulks),
    phenoData = Biobase::AnnotatedDataFrame(bulk.pheno)
  )

  if (is.null(model)) {
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
    deg.per.type <- try(
      marker_genes(
        exprs = exprs,
        pheno = pheno,
        sig.types = NULL,
        cell.type.column = cell.type.column
      )
    )

    # BSEQ-sc marker selection produces errors sometimes
    if (class(deg.per.type) == "try-error") {
      warning("Error while trying to find marker genes for BSEQ-sc")
      return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
    }
    if (length(deg.per.type) == 0) {
      warning("No genes passed the BSEQ-sc criteria")
      return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
    }

    # create single cell expression set
    sc.exprs <- Biobase::ExpressionSet(
      assayData = Matrix::as.matrix(exprs),
      phenoData = Biobase::AnnotatedDataFrame(pheno)
    )

    # create reference matrix based on single cell data and DEGs
    B <- try({
      bseqsc::bseqsc_basis(
        sc.exprs,
        deg.per.type,
        clusters = cell.type.column,
        samples = patient.column,
        ct.scale = TRUE
      )
    })
    if (length(class(B)) == 1)
    {
      if (class(B) == "try-error") {
        warning("BSEQ-sc signature matrix creation failed")
        return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
      }
    }
    
    model <- list(B = B)
  }else{
    B <- model$B
    if (!is.null(model_exclude)) {
      cts <- colnames(B)
      if (all(model_exclude %in% cts)) {
        cts <- cts[-which(cts %in% model_exclude)]
        B <- B[, cts, drop = FALSE]
      }else{
        stop("Not all cell types in 'model_exclude' are present in the model")
      }
    }
  }

  # deconvolute bulks using reference matrix B
  fit <- try(
    bseqsc::bseqsc_proportions(bulks, B, log = F, verbose = verbose)
  )
  if (length(class(fit)) == 1) {
    if (class(fit) == "try-error") {
      warning("BSEQ-sc estimation failed")
      return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
    }
  }
  

  # complete the estimation matrix in case of dropout cell types
  if (!all(unique(pheno[[cell.type.column]]) %in% rownames(fit$coefficients))) {
    est.props <- complete_estimates(
      fit$coefficients,
      unique(pheno[[cell.type.column]])
    )
  }

  return(list(
    est.props = est.props,
    sig.matrix = B,
    model = model
  ))
}
