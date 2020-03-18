# written by Tim Mirus

geneset_benchmark <- function(training.exprs, training.pheno, test.exprs, test.pheno, genesets, algorithms, bulk.data, n.repeats, exclude.from.signature = NULL){
  geneset.lists <- list()
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
      verbose = TRUE,
      split.data = FALSE,
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
