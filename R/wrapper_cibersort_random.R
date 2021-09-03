#' deconvolute given bulks with CIBERSORT using single cell data
#' (random gene selection)
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
#' @param model model for CIBERSORT deconvolution as returned by this wrapper,
#' default NULL
#' @param model_exclude character vector, cell type(s) to exclude
#' from the supplied pre-trained model, default NULL
#' @return list with four entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained\cr
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' (features x cell types)\cr
#' 3) model - list containing reference.X (signature matrix)\cr
#' @export

# source the CIBERSORT function (not available as package)
run_cibersort_random_genes <- function(
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
  suppressMessages(library(e1071, quietly = TRUE))
  suppressMessages(library(parallel, quietly = TRUE))
  suppressMessages(library(preprocessCore, quietly = TRUE))

  # error checking
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

  if (is.null(model)) {
    if (scale.cpm) {
        # prepare phenotype data and cell types to use
        exprs <- scale_to_count(exprs)
    }

    # create signature matrix
    random.genes <- sample(rownames(exprs), size = 750, replace = FALSE)
    ref.profiles <- matrix(0,
      nrow = length(random.genes),
      ncol = length(unique(pheno[[cell.type.column]]))
    )
    rownames(ref.profiles) <- random.genes
    colnames(ref.profiles) <- unique(pheno[[cell.type.column]])
    for (ct in colnames(ref.profiles)) {
        ref.profiles[,ct] <- Matrix::rowMeans(
          exprs[
            which(rownames(exprs) %in% rownames(ref.profiles)),
            which(pheno[[cell.type.column]] == ct),
            drop = F
            ]
          )
    }

    df.sig <- data.frame(GeneSymbol = rownames(ref.profiles))
    df.sig <- cbind(df.sig, ref.profiles)
    model <- list(reference.X = df.sig)
  }else{
    df.sig <- model$reference.X
    if (!is.null(model_exclude)) {
      if (all(model_exclude %in% colnames(df.sig))) {
        to_remove <- which(colnames(df.sig) %in% model_exclude)
        df.sig <- df.sig[, -to_remove, drop = FALSE]
      }else{
        stop("Not all cell types in 'model_exclude' are present in the model")
      }
    }
  }

  # CIBERSORT expects input to be supplied as .txt files
  dir.create("CIBERSORT/")
  write.table(
    df.sig, file = "CIBERSORT/signature_matrix.txt",
    quote = FALSE, row.names = FALSE, sep = "\t"
  )

  # create data frame containing bulks
  df.mix <- data.frame(GeneSymbol = rownames(bulks))
  df.mix <- cbind(df.mix, Matrix::as.matrix(bulks))
  write.table(
    df.mix, file = "CIBERSORT/mixture.txt",
    quote = FALSE, row.names = FALSE, sep = "\t"
  )

  # call CIBERSORT; quantile normalization is recommended by the authors
  # switch off permutation, as we are not interested in p-values
  result <- try({CIBERSORT(
      sig_matrix = "CIBERSORT/signature_matrix.txt",
      mixture_file = "CIBERSORT/mixture.txt",
      QN = TRUE, perm = 0
  )})
  if (class(result) == "try-error") {
    return(list(est.props = NULL, sig.matrix = NULL))
  }

  # drop the additional information in the last 3 columns
  last_cols <- (ncol(result) - 2):ncol(result)
  est.props <- t(result[seq_len(ncol(bulks)), -last_cols, drop = FALSE])

  # complete the estimation matrix in case of cell type dropouts
  if (!all(colnames(ref.profiles) %in% rownames(est.props))) {
      est.props <- complete_estimates(est.props, colnames(ref.profiles))
  }
  # CIBERSORT automatically stores the results in a file,
  # but we do not need it
  file.remove("CIBERSORT-Results.txt")
  file.remove("CIBERSORT/signature_matrix.txt")
  file.remove("CIBERSORT/mixture.txt")
  unlink("CIBERSORT", recursive = TRUE)

  return(list(est.props = est.props, sig.matrix = ref.profiles, model = model))
}
