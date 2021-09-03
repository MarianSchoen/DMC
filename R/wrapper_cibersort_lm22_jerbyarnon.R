#' deconvolute given bulks with CIBERSORT using LM22 model
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
##' @param max.genes numeric, maximum number of genes that will be included in
#' the signature for each celltype, default 500
#' @param verbose boolean
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @param scale.cpm boolean, unused
#' @return list with four entries:
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm
#' (features x cell types)
#' 3) model - NULL\cr
#' @export

# source the CIBERSORT function (not available as package)
run_cibersort_lm22_jerbyarnon <- function(
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

  if (scale.cpm) {
      # prepare phenotype data and cell types to use
      exprs <- scale_to_count(exprs)
  }

  lm22 <- read.table(
    system.file("./", "LM22.txt", package = "DAB"),
    sep = "\t", header = TRUE
  )
  lm22.genes <- as.character(lm22$Gene.symbol)
  lm22 <- as.matrix(lm22[, -1])
  rownames(lm22) <- lm22.genes

  if (length(intersect(rownames(lm22), rownames(exprs))) == 0) {
    warning(
      "No common genes in lm22 and expression matrix.
      Try converting lm22 gene names to ensembl ids..."
    )
    source("/data/tim/dab-paper-source-code/src/99_helper_functions.R")
    ensembl.genes <- convert_gene_names(lm22.genes, to = "ensembl")

    lm22 <- lm22[names(ensembl.genes), ]
    rownames(lm22) <- as.character(ensembl.genes)
  }

  lm22.genes <- intersect(rownames(lm22), rownames(exprs))
  if (length(lm22.genes) == 0) {
    warning(
      "No matching gene names in LM22 and expression matrix.
      Returning NULL."
    )
    return(est.props = NULL, sig.matrix = NULL, model = NULL)
  }

  ref.profiles <- lm22[lm22.genes, ]

  df.sig <- data.frame(GeneSymbol = rownames(ref.profiles))
  df.sig <- cbind(df.sig, ref.profiles)

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
  result <- CIBERSORT(
      sig_matrix = "CIBERSORT/signature_matrix.txt",
      mixture_file = "CIBERSORT/mixture.txt",
      QN = TRUE, perm = 0
  )

  # drop the additional information in the last 3 columns
  last_cols <- (ncol(result) - 2):ncol(result)
  est.props <- t(result[seq_len(ncol(bulks)), -last_cols, drop = FALSE])

  # complete the estimation matrix in case of cell type dropouts
  if (!all(colnames(ref.profiles) %in% rownames(est.props))) {
      est.props <- complete_estimates(est.props, colnames(ref.profiles))
  }

  # fit estimations of LM22 to available cell types
  # (hardcoded, but for now it has to suffice)
  est.props.fitted <- matrix(0, nrow = 16, ncol = ncol(est.props))
  colnames(est.props.fitted) <- colnames(est.props)
  est.props.fitted[1, ] <- colSums(est.props[c(1, 2), ]) # B
  est.props.fitted[2, ] <- colSums(est.props[c(5, 6, 7), ]) # CD4, 4 is cd8+
  est.props.fitted[3, ] <- colSums(est.props[c(11, 12), ])
  est.props.fitted[4, ] <- colSums(est.props[c(13, 14, 15, 16), ])
  used <- c(1, 2, 5, 6, 7, 11, 12, 13, 14, 15, 16)
  est.props <- est.props[-used, ]
  est.props.fitted[5:15, ] <- est.props
  rownames(est.props.fitted) <- c(
    "B.cell", "cd4+", "nk", "mono", rownames(est.props), "malignant"
  )

  # CIBERSORT automatically stores the results in a file,
  # but we do not need it
  file.remove("CIBERSORT-Results.txt")
  file.remove("CIBERSORT/signature_matrix.txt")
  file.remove("CIBERSORT/mixture.txt")
  unlink("CIBERSORT", recursive = TRUE)

  return(list(
    est.props = est.props.fitted,
    sig.matrix = ref.profiles,
    model = NULL
  ))
}
