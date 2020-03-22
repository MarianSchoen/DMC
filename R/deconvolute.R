# written by Tim Mirus

#' deconvolute given bulks with all supplied algorithms and training data, optionally plotting the results in a table plot
#'
#' @param training.expr matrix containing single-cell expression profiles (training set, one cell per column)
#' @param training.pheno data frame containing phenotype data of the single-cell training set. Has to contain column "cell_type"
#' @param test.expr matrix containing single-cell expression profiles (test set, one cell per column)
#' @param test.pheno data frame containing phenotype data of the single-cell test set. Has to contain column "cell_type"
#' @param algorithms List containing a list for each algorithm. Each sublist contains 1) name  and 2) function
#' @param plot.path character, file path/name  where the optional plot should be saved
#' @param verbose logical, default FALSE
#' @param split.data logical, if TRUE (default) then 10% of the training data will be used for reference profile creation and
#' the rest for feature selection/optimization
#' @param exclude.from.bulks character vector containing cell types to be excluded from the bulks (if they are not supplied).
#' If not specified, all will be used.
#' @param exclude.from.signature character vector containing cell types to be excluded from the signature matrix.
#' If not specified, all will be used.
#' @param optimize logical, should the signature matrix be optimized by condition number? If FALSE, max.genes genes will be used
#' @param max.genes maximum number of genes that will be included in the signature for each celltype
#' @param cor.title character, title of the optional table plot
#' @param time.title character, title of the optional runtime plot
#' @param n.bulks number of bulks to build if they are not supplied to the function, default 50
#' @param plot.results logical, should performance and runtime plots be created? default TRUE
#' @param bulks matrix containing expression profiles of bulks in the columns. If not supplied, bulks will be created
#' @param scatter.plots boolean, if TRUE create scatter plots of real vs estimated quantities of cell types for all algorithms; default FALSE
#' @param n.repeats integer determining the number of times deconvolution should be repeated for each algorithm, default 1
#' @return list with two entries:
#' 1) est.props: matrix containing estimated proportions / quantities for all cell types in all bulks (cell type x bulk)
#' 2) bulk.props: matrix containing the real proportions / quantities for all cell types in all bulks (cell type x bulk)
#' @example deconvolute(training.expr, training.pheno, test.expr, test.pheno, algorithm_list, "results/test", cor.title = 'correlation plot', n.repeats = 10)

suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

deconvolute <- function(training.expr,
                        training.pheno,
                        test.expr,
                        test.pheno,
                        algorithms,
                        verbose = FALSE,
                        split.data = FALSE,
                        exclude.from.bulks = NULL,
                        exclude.from.signature = NULL,
                        optimize = TRUE,
                        max.genes = 500,
                        n.bulks = 500,
                        bulks = NULL,
                        n.repeats = 1) {
  # parameter checks
  if (n.repeats < 1) {
    n.repeats <- 1
  }
  if (nrow(training.pheno) != ncol(training.expr)) {
      stop("Number of columns in training.expr and rows in training.pheno do not match")
  }
  if(!is.null(test.pheno) && !is.null(test.exprs)){
  	if (nrow(test.pheno) != ncol(test.expr)) {
      		stop("Number of columns in test.expr and rows in test.pheno do not match")
  	}
  }else{
	  if(is.null(bulks)){
		  stop("Have to supply either test data or pre-built bulks")
	  }
  }
  if (!is.null(max.genes) && max.genes == 0) {
      max.genes <- NULL
  }
  if (!length(algorithms) > 0 || any(sapply(algorithms, length) != 2)) {
    stop("Check algorithm list")
  }
  if (is.null(bulks)) {
    # create validation bulks with known proportions
    if (verbose)
      print("Creating bulks")
    # remove cells whose type should not be in the bulks
    if (!is.null(exclude.from.bulks)) {
      if (length(which(unique(test.pheno[, "cell_type"]) %in% exclude.from.bulks)) > 0) {
        include.in.bulks <-
          unique(test.pheno[, "cell_type"])[-which(unique(test.pheno[, "cell_type"]) %in% exclude.from.bulks)]
      }else{
        include.in.bulks <- unique(test.pheno[, "cell_type"])
      }
    } else {
      include.in.bulks <- unique(test.pheno[, "cell_type"])
    }
    # if seeding is set the bulks will always be the same
    if (seeding)
      set.seed(10)
    validation.data <- create_bulks(test.expr,
                                    test.pheno,
                                    n.bulks,
                                    include.in.bulks,
                                    sum.to.count = TRUE)
    real.props <- validation.data$props
    bulks.expr <- validation.data$bulks
  } else {
    # use the supplied bulks
    real.props <- bulks$props
    bulks.expr <- bulks$bulks

    if (nrow(training.expr) != nrow(bulks.expr)) {
      features <- intersect(rownames(training.expr), rownames(bulks.expr))
      if (length(features) > 0) {
          training.expr <- training.expr[features, ]
          bulks.expr <- bulks.expr[features, ]
      }
  }
  }
  
  # deconvolute bulks with all supplied algorithms
  if (verbose)
    print("Deconvoluting...")
  results.list <- list()
  # perform deconvolution several times
  for (i in seq_len(n.repeats)) {
    if (verbose && n.repeats > 1)
      cat("Repetition ", i, "\n", sep = "")
    results.list[[as.character(i)]] <- list()
    for (f in algorithms) {
      if (verbose)
        print(f$name)
      time <- system.time({
        results.list[[as.character(i)]][[f$name]] <- f$algorithm(
          training.expr,
          training.pheno,
          bulks.expr,
          exclude.from.signature,
          max.genes = max.genes,
          optimize = optimize,
          split.data = split.data
        )
      })[3]
      results.list[[as.character(i)]][[f$name]]$name <- f$name
      results.list[[as.character(i)]][[f$name]]$times <- time
    }
  }
  return(list(results.list = results.list, bulk.props = real.props))
  }
