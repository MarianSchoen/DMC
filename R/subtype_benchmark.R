#' perform deconvolution benchmark on simulated bulks with simulated subtypes
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
#' @param verbose logical, default FALSE
#' @param split.data logical, if TRUE (default) then 10% of the training data will be used for reference profile creation and
#' the rest for feature selection/optimization
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? 
#' @return list containing deconvolution results

subtype_benchmark <- function(training.exprs, 
                              training.pheno, 
                              test.exprs, 
                              test.pheno, 
                              algorithms, 
                              bulk.data, 
                              n.repeats, 
                              exclude.from.signature = NULL, 
                              verbose = F, 
                              split.data = FALSE,
                              cell.type.column = "cell_type"){
  # parameter checks
  if(ncol(training.exprs) != nrow(training.pheno)){
    stop("training.exprs and training.pheno do not match")
  }
  if(!is.null(test.exprs) || !is.null(test.pheno)){
    stop("test.exprs and test.pheno are not being used in this benchmark and must be NULL")
  }
  if(!is.null(bulk.data)){
    if(!is.list(bulk.data) || !c("bulks", "props") %in% names(bulk.data)){
      stop("bulk.data has the wrong format")
    }
  }
  if(!"subtype" %in% colnames(training.pheno) || !cell.type.column %in% colnames(training.pheno)){
    stop("Required columns missing from training.pheno")
  }
  if(!is.numeric(n.repeats)){
    stop("n.repeats has to be numeric")
  }
  if(n.repeats < 1){
    warning("n.repeats has to be greater than 0. setting to 1")
    n.repeats <- 1
  }
  # algorithm list is checked in deconvolute()

  # replace cell_type with combination of cell_type and subtype,
  # then add original cell types as column coarse_type to the pheno data
  # call deconvolute with 'subtype=TRUE'
  temp.pheno <- training.pheno
  temp.pheno[, cell.type.column] <- paste(temp.pheno[,cell.type.column], temp.pheno[,"subtype"], sep = ".")
  temp.pheno <- cbind(temp.pheno, coarse_type = training.pheno[[cell.type.column]])
  result <- deconvolute(training.exprs, temp.pheno, NULL, NULL, algorithms, verbose, split.data, NULL, exclude.from.signature, bulks = list(bulks = bulk.data$bulks, props = bulk.data$props), n.repeats = n.repeats, subtypes = TRUE, cell.type.column = cell.type.column)
  return(result)
}
