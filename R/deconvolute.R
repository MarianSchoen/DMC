#' deconvolute given bulks with all supplied algorithms and training data
#'
#' @param training.expr matrix containing single-cell expression profiles
#' (training set, one cell per column)
#' @param training.pheno data frame containing phenotype data of the single-cell
#'  training set. Has to contain column 'cell.type.column'
#' @param test.expr matrix containing single-cell expression profiles
#' (test set, one cell per column)
#' @param test.pheno data frame containing phenotype data of the single-cell
#'  test set. Has to contain column 'cell.type.column'
#' @param algorithms List containing a list for each algorithm.
#' Each sublist contains 1) name \cr 2) function \cr 3) model
#' @param verbose logical, default FALSE
#' @param exclude.from.bulks character vector containing cell types to be
#' excluded from the bulks (if they are not supplied).
#' If not specified, all will be used.
#' @param exclude.from.signature character vector containing cell types to be
#' excluded from the signature matrix.
#' If not specified, all will be used.
#' @param max.genes maximum number of genes that will be included in the
#' signature for each celltype
#' @param n.bulks number of bulks to build if they are not supplied to the
#' function, default 500
#' @param bulks matrix containing expression profiles of bulks in the columns.
#' If not supplied, bulks will be created
#' @param n.repeats integer determining the number of times deconvolution
#' should be repeated for each algorithm, default 1
#' @param cell.type.column string, which column of 'training.pheno'/'test.pheno'
#' holds the cell type information? default 'cell_type'
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default 'patient'
#' @param subtypes boolean, are simulated subtypes used for deconvolution?
#' @param n.profiles.per.bulk positive numeric, number of samples to be randomly
#' drawn for each simulated bulk; default 1000; only needed when bulks=NULL
#'
#' @return list with two entries:
#' 1) results.list: list containing deconvolution results for all algorithms
#' and repetitions as returned by the algorithm functions
#' 2) bulk.props: matrix containing the real proportions / quantities for all
#' cell types in all bulks (cell type x bulk)
#' @export
deconvolute <- function(
  training.expr,
  training.pheno,
  test.expr,
  test.pheno,
  algorithms,
  verbose = FALSE,
  exclude.from.bulks = NULL,
  exclude.from.signature = NULL,
  max.genes = 500,
  n.bulks = 500,
  bulks = NULL,
  n.repeats = 1,
  subtypes = FALSE,
  cell.type.column = "cell_type",
  patient.column = "patient",
  n.profiles.per.bulk = 1000
  ) {
  # parameter checks
  if (n.repeats < 1) {
    n.repeats <- 1
  }
  if (nrow(training.pheno) != ncol(training.expr) && !is.null(training.expr)) {
      stop("Number of columns in training.expr and
           rows in training.pheno do not match")
  }
  if (!is.null(test.pheno) && !is.null(test.expr)) {
  	if (nrow(test.pheno) != ncol(test.expr)) {
      		stop("Number of columns in test.expr and
               rows in test.pheno do not match")
  	}
  }else{
	  if (is.null(bulks)) {
		  stop("Need either test data or pre-built bulks")
	  }
  }
  if (!is.null(max.genes) && max.genes == 0) {
      max.genes <- NULL
  }

  if (!length(algorithms) > 0 || any(sapply(algorithms, length) < 3)) {
    stop("Check algorithm list")
  }
  if (is.null(bulks)) {
    # create validation bulks with known proportions
    if (verbose) {
      cat("creating artificial bulks for simulation\n")
    }
    # remove cells whose type should not be in the bulks
    cts <- unique(test.pheno[, cell.type.column])
    include.in.bulks <- cts
    if (!is.null(exclude.from.bulks)) {
      if (length(which(cts %in% exclude.from.bulks)) > 0) {
        include.in.bulks <- cts[-which(cts %in% exclude.from.bulks)]
      }
    }
    set.seed(1)
    validation.data <- create_bulks(
      exprs = test.expr
      , pheno = test.pheno
      , cell.type.column = cell.type.column
      , n.bulks = n.bulks
      , include.in.bulks = include.in.bulks
      , sum.to.count = TRUE
      , n.profiles.per.bulk = n.profiles.per.bulk
      , frac = 10
    )
    real.props <- validation.data$props
    bulks.expr <- validation.data$bulks
  }else{
    # use the supplied bulks
    real.props <- bulks$props
    bulks.expr <- bulks$bulks
    if (!is.null(training.expr)) {
      if (nrow(training.expr) != nrow(bulks.expr)) {
        features <- intersect(rownames(training.expr), rownames(bulks.expr))
        if (length(features) > 0) {
            training.expr <- training.expr[features, ]
            bulks.expr <- bulks.expr[features, ]
        }
      }
    }
  }

  # deconvolute bulks with all supplied algorithms
  if (verbose) {
    cat("deconvolving:\n")
  }
  results.list <- list()
  # perform deconvolution several times
  for (i in seq_len(n.repeats)) {
    if (verbose && n.repeats > 1) {
      cat("Repetition ", i, " of ", n.repeats, "\n", sep = "")
    }
    results.list[[as.character(i)]] <- list()
    for (f in algorithms) {
      if (verbose) {
        cat(f$name, "\t", sep = "")
        cat(i, "\n")
      }
      time <- system.time({
        results.list[[as.character(i)]][[f$name]] <- f$algorithm(
          training.expr,
          training.pheno,
          bulks.expr,
          exclude.from.signature,
          cell.type.column = cell.type.column,
          max.genes = max.genes,
          patient.column = patient.column
        )
      })[3]
      results.list[[as.character(i)]][[f$name]]$name <- f$name
      results.list[[as.character(i)]][[f$name]]$times <- time
    }
    if (verbose) cat("\n")
  }

  return(list(results.list = results.list, bulk.props = real.props))
}
