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
#' @param algorithms List containing a list for each algorithm. Each sublist contains 1) name  and 2) function
#' @param verbose logical, default FALSE
#' @param split.data logical, if TRUE (default) then 10 \% of the training data will be used for reference profile creation and
#' the rest for feature selection/optimization 
#' @param exclude.from.bulks character vector containing cell types to be excluded from the bulks (if they are not supplied).
#' If not specified, all will be used.
#' @param exclude.from.signature character vector containing cell types to be excluded from the signature matrix.
#' If not specified, all will be used.
#' @param optimize logical, should the signature matrix be optimized by condition number? If FALSE, max.genes genes will be used
#' @param max.genes maximum number of genes that will be included in the signature for each celltype
#' @param n.bulks number of bulks to build if they are not supplied to the function, default 500
#' @param bulks matrix containing expression profiles of bulks in the columns. If not supplied, bulks will be created
#' @param n.repeats integer determining the number of times deconvolution should be repeated for each algorithm, default 1
#' @param cell.type.column string, which column of 'training.pheno'/'test.pheno'
#' holds the cell type information? 
#' @param subtypes boolean, are simulated subtypes used for deconvolution?
#'
#' @return list with two entries:
#' 1) results.list: list containing deconvolution results for all algorithms and repetitions as returned by the algorithm functions
#' 2) bulk.props: matrix containing the real proportions / quantities for all cell types in all bulks (cell type x bulk)

deconvolute <- function(
  training.expr,
  training.pheno,
  test.expr,
  test.pheno,
  cell.type.column = "cell_type", 
  algorithms,
  verbose = FALSE,
  split.data = FALSE,
  exclude.from.bulks = NULL,
  exclude.from.signature = NULL,
  optimize = TRUE,
  max.genes = 500,
  n.bulks = 500,
  bulks = NULL,
  n.repeats = 1,
  subtypes = FALSE
  ) {
  # parameter checks
  if (n.repeats < 1) {
    n.repeats <- 1
  }
  if (nrow(training.pheno) != ncol(training.expr)) {
      stop("Number of columns in training.expr and rows in training.pheno do not match")
  }
  if(!is.null(test.pheno) && !is.null(test.expr)){
  	if (nrow(test.pheno) != ncol(test.expr)) {
      		stop("Number of columns in test.expr and rows in test.pheno do not match")
  	}
  }else{
	  if(is.null(bulks)){
		  stop("Need either test data or pre-built bulks")
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
      cat("creating artificial bulks for simulation\n")
    # remove cells whose type should not be in the bulks
    if (!is.null(exclude.from.bulks)) {
      if (length(which(unique(test.pheno[, cell.type.column]) %in% exclude.from.bulks)) > 0) {
        include.in.bulks <-
          unique(test.pheno[, cell.type.column])[-which(unique(test.pheno[, cell.type.column]) %in% exclude.from.bulks)]
      }else{
        include.in.bulks <- unique(test.pheno[, cell.type.column])
      }
    } else {
      include.in.bulks <- unique(test.pheno[, cell.type.column])
    }
    validation.data <- create_bulks(
      exprs = test.expr
      , pheno = test.pheno
      , cell.type.column = cell.type.column 
      , n.bulks = n.bulks
      , include.in.bulks = include.in.bulks
      , sum.to.count = TRUE
      )
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
    cat("deconvolving:\n")
  results.list <- list()
  # perform deconvolution several times
  for (i in seq_len(n.repeats)) {
    if (verbose && n.repeats > 1)
      cat("Repetition ", i, " of ", n.repeats, "\n", sep = "")
    results.list[[as.character(i)]] <- list()
    for (f in algorithms) {
      if (verbose)
        cat(f$name, "\t", sep = "")
      time <- system.time({
        results.list[[as.character(i)]][[f$name]] <- f$algorithm(
          training.expr,
          training.pheno,
          bulks.expr,
          exclude.from.signature,
          cell.type.column = cell.type.column,
          max.genes = max.genes,
          optimize = optimize,
          split.data = split.data
        )
      })[3]
      results.list[[as.character(i)]][[f$name]]$name <- f$name
      results.list[[as.character(i)]][[f$name]]$times <- time
      
      if(subtypes){
        if(! "coarse_type" %in% colnames(training.pheno)){
          stop("'subtypes' is TRUE, but column 'coarse_type' is missing from pheno data")
        }
        # change the est.props to coarse types by combining subtypes
        if(!is.null(results.list[[as.character(i)]][[f$name]]$est.props)){
        	coarse.rnames <- sapply(strsplit(rownames(results.list[[as.character(i)]][[f$name]]$est.props), ".", fixed = TRUE), function(x){x[1]})
        	temp.props <- matrix(0, nrow = length(unique(coarse.rnames)), ncol = ncol(results.list[[as.character(i)]][[f$name]]$est.props))
        	rownames(temp.props) <- unique(coarse.rnames)
        	colnames(temp.props) <- colnames(bulks.expr)
        	for(t in rownames(temp.props)){
          	idx <- which(coarse.rnames == t)
          	temp.props[t,] <- colSums(results.list[[as.character(i)]][[f$name]]$est.props[idx,,drop=F])
        	}
        	results.list[[as.character(i)]][[f$name]]$est.props <- temp.props
      	}
      }
    }
    if(verbose) cat("\n")
  }

  return(list(results.list = results.list, bulk.props = real.props))
  }
