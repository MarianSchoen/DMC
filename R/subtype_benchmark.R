subtype_benchmark <- function(training.exprs, training.pheno, test.exprs, test.pheno, algorithms, bulk.data, n.repeats, exclude.from.signature = NULL, verbose = F){
  temp.pheno <- training.pheno
  temp.pheno[, "cell_type"] <- paste(temp.pheno[,"cell_type"], temp.pheno[,"subtype"], sep = ".")
  temp.pheno <- cbind(temp.pheno, coarse_type = training.pheno$cell_type)
  result <- deconvolute(training.exprs, temp.pheno, NULL, NULL, algorithms, verbose, FALSE, NULL, exclude.from.signature, bulks = list(bulks = bulk.data$bulks, props = bulk.data$props), n.repeats = n.repeats, subtypes = TRUE)
}