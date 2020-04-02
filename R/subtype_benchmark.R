subtype_benchmark <- function(training.exprs, training.pheno, test.exprs, test.pheno, algorithms, bulk.data, n.repeats, exclude.from.signature = NULL, verbose = F, split.data = FALSE){
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
  if(!"subtype" %in% colnames(training.pheno) || !"cell_type" %in% colnames(training.pheno)){
    stop("Required columns missing from training.pheno")
  }
  if(!is.numeric(n.repeats)){
    stop("n.repeats has to be numeric")
  }
  if(n.repeats < 1){
    warning("n.repeats has to be greater than 0. setting to 1")
    n.repeats <- 1
  }
  # algorithm list will be checked in deconvolute()

  temp.pheno <- training.pheno
  temp.pheno[, "cell_type"] <- paste(temp.pheno[,"cell_type"], temp.pheno[,"subtype"], sep = ".")
  temp.pheno <- cbind(temp.pheno, coarse_type = training.pheno$cell_type)
  result <- deconvolute(training.exprs, temp.pheno, NULL, NULL, algorithms, verbose, split.data, NULL, exclude.from.signature, bulks = list(bulks = bulk.data$bulks, props = bulk.data$props), n.repeats = n.repeats, subtypes = TRUE)
}
