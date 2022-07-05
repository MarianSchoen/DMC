#' perform deconvolution of simulated bulks on different gene sets
#'
#' @param training.exprs matrix containing single-cell expression profiles
#' (training set, one cell per column)
#' @param training.pheno data frame containing phenotype data of the
#' single-cell training set. Has to contain column `cell.type.column`
#' @param test.exprs matrix containing single-cell expression profiles
#' (test set, one cell per column)
#' @param test.pheno data frame containing phenotype data of the
#' single-cell test set. Has to contain column `cell.type.column`
#' @param genesets list of gene sets (character vectors)
#' @param algorithms List containing a list for each algorithm.
#' Each sublist contains 1) name, 2) function and 3) model
#' @param verbose logical, default FALSE
#' @param exclude.from.signature character vector containing cell types to be
#' excluded from the signature matrix. If not specified, all will be used.
#' @param bulk.data list with two entries:\cr
#' 1) bulks - matrix containing expression data of the bulks
#' (one bulk per column)\cr
#' 2) props - matrix containing the true fractions of cell types within the
#' bulks (cell type x bulk)
#' @param n.repeats integer determining the number of times deconvolution
#' should be repeated for each algorithm, default 3
#' @param cell.type.column string, which column of 'training.pheno'/'test.pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default "patient"
#' @return list containing deconvolution results for all algorithms
#' for all genesets

geneset_benchmark <- function(
  training.exprs,
  training.pheno,
  test.exprs,
  test.pheno,
  genesets,
  algorithms,
  bulk.data,
  n.repeats = 3,
  exclude.from.signature = NULL,
  verbose = FALSE,
  cell.type.column = "cell_type",
  patient.column = "patient"
) {
  # parameter checks
  gene_intersects <- sapply(genesets, function(x) {
    any(x %in% rownames(training.exprs))
  })
  if (!all(gene_intersects)) {
    stop("one or more genesets do not contain
          any genes present in the expression data")
  }
  if (ncol(training.exprs) != nrow(training.pheno)) {
    stop("training.exprs and training.pheno do not match")
  }
  if (!is.null(test.exprs) || !is.null(test.pheno)) {
    if (ncol(test.exprs) != nrow(test.pheno)) {
      stop("test.exprs and test.pheno do not match")
    }
  }
  if (!is.null(bulk.data)) {
    if (!is.list(bulk.data) || !c("bulks", "props") %in% names(bulk.data)) {
      stop("bulk.data has the wrong format")
    }
  }

  geneset.lists <- list()
  # deconvolve using each geneset
  for (i in seq_len(length(genesets))) {
    genes <- genesets[[i]]
    geneSetName <- names(genesets)[i]
    if (!is.na(as.numeric(geneSetName))) {
	    geneSetName <- paste0("geneset_", geneSetName)
    }
   
    # reduce to genes contained in current gene set
    temp.exprs <- training.exprs[which(rownames(training.exprs) %in% genes), ]

    # deconvolve n.repeats times for each gene set
    temp.results <- deconvolute(
      training.expr = temp.exprs,
      training.pheno = training.pheno,
      test.expr = test.exprs,
      test.pheno = test.pheno,
      algorithms = algorithms,
      verbose = verbose,
      exclude.from.signature = exclude.from.signature,
      max.genes = nrow(temp.exprs),
      n.bulks = 0,
      bulks = bulk.data,
      n.repeats = n.repeats,
      cell.type.column = cell.type.column,
      patient.column = patient.column
    )
    # add only the results, not the real proportions that are returned
    geneset.lists[[geneSetName]] <- temp.results[[1]]
  }
  # add real props once at the top level
  geneset.lists[["bulk.props"]] <- temp.results$bulk.props
  return(geneset.lists)
}
