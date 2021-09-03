#' training set size benchmark: deconvolute given bulks using training
#' sets of different sizes to train deconvolution models
#'
#' @param training.exprs matrix containing single-cell expression profiles
#' (training set, one cell per column)
#' @param training.pheno data frame containing phenotype data of the
#' single-cell training set. Has to contain column "cell_type"
#' @param test.exprs matrix containing single-cell expression profiles
#' (test set, one cell per column)
#' @param test.pheno data frame containing phenotype data of the single-cell
#' test set. Has to contain column `cell.type.column`
#' @param algorithms List containing a list for each algorithm.
#' Each sublist contains 1) name, 2) function and 3) model
#' @param bulk.data list with two entries:\cr
#' 1) bulks - matrix containing expression data of the bulks
#' (one bulk per column)\cr
#' 2) props - matrix containing the true fractions of cell
#' types within the bulks (cell type x bulk)
#' @param n.repeats integer determining the number of times deconvolution
#' should be repeated for each algorithm
#' @param exclude.from.signature character vector containing cell types to be
#' excluded from the signature matrix. If not specified, all will be used.
#' @param step.size numerical 0 < step.size < 1; fraction of samples
#'  by which size of training set is increased each step; default 0.05
#' @param verbose logical, default FALSE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default "patient"
#' @return list containing deconvolution results for all algorithms
#' for each training set size

sample_size_benchmark <- function(
  training.exprs,
  training.pheno,
  test.exprs,
  test.pheno,
  algorithms,
  bulk.data,
  n.repeats,
  exclude.from.signature = NULL,
  step.size = 0.05,
  verbose = FALSE,
  cell.type.column = "cell_type",
  patient.column = "patient"
) {
  # parameter checks
  if (ncol(training.exprs) != nrow(training.pheno)) {
      stop("training.exprs and training.pheno do not match")
  }
  if (!is.null(test.exprs) || !is.null(test.pheno)) {
    if (ncol(test.exprs) != nrow(test.pheno)) {
      stop("test.exprs and test.pheno do not match")
    }
  }
  if (!is.null(bulk.data)) {
    if (!is.list(bulk.data) || !all(c("bulks", "props") %in% names(bulk.data))) {
      stop("bulk.data has the wrong format")
    }
  }
  if (!is.numeric(n.repeats)) {
    stop("n.repeats has to be numeric")
  }
  if (n.repeats < 1) {
    warning("n.repeats has to be greater than 0. setting to 1")
    n.repeats <- 1
  }
  if (!is.numeric(step.size) || step.size <= 0 || step.size >= 1) {
    stop("step.size must be numeric in (0,1)")
  }

  # list for storing the results for each fraction
  sample.size.lists <- list()
  for (i in seq(1, 1 / step.size)) {
    sample.size.lists[[as.character(i * step.size)]] <- list()
  }
  cell.types <- unique(training.pheno[, cell.type.column])

  # repeats n.repeats times
  for (rep in seq_len(n.repeats)) {
    if (verbose) cat("Repetition ", rep, " of ", n.repeats, "\n", sep = "")
    # a pool of available samples
    available.samples <- seq_len(ncol(training.exprs))

    # growing sc data structures
    temp.expr <- c()
    temp.pheno <- c()

    # deconvolve with growing training set
    for (step in seq(1, 1 / step.size)) {
      if (verbose) cat(step * step.size * 100, "%\t", sep = "")
      # grow training set by certain amount (randomly selected) for each type
      for (t in cell.types) {
        samples.to.add <- c()
        # do not try to sample more cells than there are left of this type
        n.cells <- min(ceiling(
          step.size * length(which(training.pheno[, cell.type.column] == t))),
          length(which(training.pheno[available.samples, cell.type.column] == t)
        ))
        # check if there are any unused cells of this type
        # sample n.cells samples for this cell type
        if (any(training.pheno[available.samples, cell.type.column] == t)) {
          samples.to.add <- sample(
            which(training.pheno[available.samples, cell.type.column] == t),
            size = n.cells,
            replace = FALSE
          )
        }
        # extend the set
        if (length(samples.to.add) > 0) {
          temp.expr <- cbind(
            temp.expr,
            training.exprs[, available.samples[samples.to.add], drop = FALSE]
          )
          temp.pheno <- rbind(
            temp.pheno,
            training.pheno[available.samples[samples.to.add], ]
          )
          # remove the used samples from the pool
          available.samples <- available.samples[-samples.to.add]
        }
      }

      # deconvolve once with this training set
      temp.results <- deconvolute(
        training.expr = temp.expr,
        training.pheno = temp.pheno,
        test.expr = test.exprs,
        test.pheno = test.pheno,
        algorithms = algorithms,
        verbose = FALSE,
        exclude.from.signature = exclude.from.signature,
        max.genes = NULL,
        n.bulks = 0,
        bulks = bulk.data,
        n.repeats = 1,
        cell.type.column = cell.type.column,
        patient.column = patient.column
      )
      # add only the result list returned by deconvolute() to sample.size.list
      idx <- as.character(step * step.size)
      sample.size.lists[[idx]][[as.character(rep)]] <- temp.results[[1]][[1]]
    }
    if (verbose) cat("\n")
  }
  # add bulk props to the top level of the result list
  sample.size.lists[["bulk.props"]] <- temp.results$bulk.props
  return(sample.size.lists)
}
