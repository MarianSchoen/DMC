geneset_benchmark <- function(training.exprs, training.pheno, test.exprs, test.pheno, genesets, algorithms, bulk.data, n.repeats, exclude.from.signature = NULL, verbose = FALSE, split.data = FALSE){
  geneset.lists <- list()
  if(!all(sapply(genesets, function(x){any(x %in% rownames(training.exprs))}))){
    stop("one or more genesets do not contain any genes present in the expression data")
  }
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
  # deconvolute using each geneset
  for (i in 1:length(genesets)) {
    genes <- genesets[[i]]
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
      split.data = split.data,
      exclude.from.signature = exclude.from.signature,
      optimize = TRUE,
      max.genes = nrow(temp.exprs),
      n.bulks = 0,
      bulks = bulk.data,
      n.repeats = n.repeats
    )
    # add only the results, not the real proportions that are returned
    geneset.lists[[names(genesets)[i]]] <- temp.results[[1]]
  }
  # add real props once at the top level
  geneset.lists[["bulk.props"]] <- temp.results$bulk.props
  return(geneset.lists)
}
