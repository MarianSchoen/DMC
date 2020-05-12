#' training set size benchmark
#' 
#' @param training.exprs matrix containing single-cell expression profiles (training set, one cell per column)
#' @param training.pheno data frame containing phenotype data of the single-cell training set. Has to contain column "cell_type"
#' @param test.exprs matrix containing single-cell expression profiles (test set, one cell per column)
#' @param test.pheno data frame containing phenotype data of the single-cell test set. Has to contain column "cell_type"
#' @param algorithms List containing a list for each algorithm. Each sublist contains 1) name  and 2) function
#' @param bulk.data list with two entries:
#' 1) bulks - matrix containing expression data of the bulks (one bulk per column)
#' 2) props - matrix containing the true fractions of cell types within the bulks (cell type x bulk)
#' @param n.repeats integer determining the number of times deconvolution should be repeated for each algorithm
#' @param exclude.from.signature character vector containing cell types to be excluded from the signature matrix.
#' If not specified, all will be used.
#' @param step.size numerical 0 < step.size < 1; fraction of samples by which size of training set is increased each step
#' @param verbose logical, default FALSE
#' @param split.data logical, if TRUE (default) then 10% of the training data will be used for reference profile creation and
#' the rest for feature selection/optimization
#' @return list containing deconvolution results for all algorithms for each training set size

sample_size_benchmark <- function(training.exprs, training.pheno, test.exprs, test.pheno, algorithms, bulk.data, n.repeats, exclude.from.signature = NULL, step.size = 0.05, verbose = FALSE, split.data = FALSE) {
# parameter checks
if(ncol(training.exprs) != nrow(training.pheno)){
    stop("training.exprs and training.pheno do not match")
  }
  if(!is.null(test.exprs) || !is.null(test.pheno)){
    if(ncol(test.exprs) != nrow(test.pheno)){
      stop("test.exprs and test.pheno do not match")
    }
  }
  if(!is.null(bulk.data)){
    if(!is.list(bulk.data) || !c("bulks", "props") %in% names(bulk.data)){
      stop("bulk.data has the wrong format")
    }
  }
  if(!is.numeric(n.repeats)){
    stop("n.repeats has to be numeric")
  }
  if(n.repeats < 1){
    warning("n.repeats has to be greater than 0. setting to 1")
    n.repeats <- 1
  }
  if(!is.numeric(step.size) || step.size <= 0 || step.size >= 1){
    stop("step.size must be numeric in (0,1)")
  }

# list for storing the results for each fraction
sample.size.lists <- list()
for(i in seq(1, 1 / step.size)){
  sample.size.lists[[as.character(i*step.size)]] <- list()
}

cell.types <- unique(training.pheno[, "cell_type"])

# repeats n.repeats times
for (rep in seq_len(n.repeats)) {
  # a pool of available samples
  available.samples <- 1:ncol(training.exprs)

  # growing sc data structures
  temp.expr <- c()
  temp.pheno <- c()

  # deconvolve with growing training set
  for (step in seq(1, 1 / step.size)) {
    # grow training set by certain amount (randomly selected) for each type
    for (t in cell.types) {
      samples.to.add <- c()
      # do not try to sample more cells than there are left of this type
      n.cells <- min(
        ceiling(step.size * length(which(training.pheno[, "cell_type"] == t))),
        length(which(training.pheno[available.samples, "cell_type"] == t))
      )
      # check if there are any unused cells of this type
      # sample n.cells samples for this cell type
      if (any(training.pheno[available.samples, "cell_type"] == t)) {
          samples.to.add <- sample(
            which(training.pheno[available.samples, "cell_type"] == t),
            size = n.cells,
            replace = FALSE
      )
      }
      # extend the set
      if (length(samples.to.add) > 0) {
        temp.expr <- cbind(
          temp.expr,
          training.exprs[, available.samples[samples.to.add], drop = F]
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
      verbose = verbose,
      split.data = split.data,
      exclude.from.signature = exclude.from.signature,
      optimize = TRUE,
      max.genes = NULL,
      n.bulks = 0,
      bulks = bulk.data,
      n.repeats = 1
    )
    # add only the result list returned by deconvolute() to sample.size.list
    sample.size.lists[[as.character(step*step.size)]][[as.character(rep)]] <- temp.results[[1]][[1]]
  }
}
# add bulk props to the top level of the result list
sample.size.lists[["bulk.props"]] <- temp.results$bulk.props
return(sample.size.lists)
}

