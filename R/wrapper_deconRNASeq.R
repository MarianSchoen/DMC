#' deconvolute given bulks with DeconRNASeq using single cell data
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
#' holds the patient information; optional, default NULL
#' @param scale.cpm boolean, scale single-cell profiles to CPM? default FALSE
#' @param model model for DeconRNASeq deconvolution as returned by this wrapper,
#' default NULL
#' @param model_exclude character vector, cell type(s) to exclude
#' from the supplied pre-trained model, default NULL
#' @return list with two entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained\cr
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' (features x cell types)\cr
#' 3) model - list containing reference.X (signature matrix)\cr
#' @example run_deconrnaseq(training.exprs, training.pheno, bulk.exprs)
#' @export
run_deconrnaseq <- function(
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
  if (is.null(model)) {
    if (is.null(exprs) || is.null(pheno)){
      stop("If no model is given, expression and pheno data are required.")
    }
    if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
    }
    features <- intersect(rownames(exprs), rownames(bulks))
    if (length(features) > 0) {
      exprs <- exprs[features, ]
      bulks <- bulks[features, ]
    }
    if (!is.null(max.genes) && max.genes == 0) {
      max.genes <- NULL
    }
    if (scale.cpm) {
      # prepare phenotype data and cell types to use
      exprs <- scale_to_count(exprs)
    }

    # create signature matrix (DeconRNASeq needs data frames)
    ref.profiles <- create_sig_matrix(exprs,
        pheno,
        exclude.from.signature,
        max.genes = max.genes,
        cell.type.column = cell.type.column
    )
    if (is.null(ref.profiles)) {
      return(list(est.props = NULL, sig.matrix = NULL, model = NULL))
    }
    model <- list(reference.X = ref.profiles)
  }else{
    ref.profiles <- model$reference.X
    if (!is.null(model_exlucde)) {
      cts <- colnames(ref.profiles)
      if (all(model_exclude %in% cts)) {
        cts <- cts[-which(cts %in% model_exclude)]
        ref.profiles <- ref.profiles[, cts, drop = FALSE]
      }else{
        stop("Not all cell types in 'model_exclude' are present in the model")
      }
    }
  }

  # create bulk data frame
  df.mix <- as.data.frame(Matrix::as.matrix(bulks))
  rownames(df.mix) <- rownames(bulks)

  # there is no option to switch the output of this function off
  # deconvolute
  invisible(capture.output(
    result <- try(
      DeconRNASeq::DeconRNASeq(df.mix, as.data.frame(ref.profiles)),
      silent = TRUE
    )
  ))

  if (length(class(result)) == 1)
  {
    if (class(result) == "try-error") {
      return(list(est.props = NULL, sig.matrix = ref.profiles, model = model))
    }
  }
  # select the interesting rows and transpose to have bulks = columns
  result <- t(result$out.all[1:ncol(bulks), , drop = FALSE])
  colnames(result) <- colnames(bulks)

  # complete estimation matrix in case of droput cell types
  if (!all(colnames(ref.profiles) %in% rownames(result))) {
    result <- complete_estimates(result, colnames(ref.profiles))
  }
  return(list(est.props = result, sig.matrix = ref.profiles, model = model))
}
