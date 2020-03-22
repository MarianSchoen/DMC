# written by Tim Mirus

sample_size_benchmark <- function(training.exprs, training.pheno, test.exprs, test.pheno, algorithms, bulk.data, n.repeats, exclude.from.signature = NULL, step.size = 0.05, verbose = FALSE) {
sample.size.lists <- list()
for(i in seq(1, 1 / step.size)){
  sample.size.lists[[as.character(i*step.size)]] <- list()
}
cell.types <- unique(training.pheno[, "cell_type"])
for (rep in seq_len(n.repeats)) {
  available.samples <- 1:ncol(training.exprs)
  temp.expr <- c()
  temp.pheno <- c()
  # deconvolve with growing training set
  for (step in seq(1, 1 / step.size)) {
    # grow training set
    for (t in cell.types) {
      samples.to.add <- c()
      # do not try to sample more cells than there are left of this type
      n.cells <- min(
        ceiling(step.size * length(which(training.pheno[, "cell_type"] == t))),
        length(which(training.pheno[available.samples, "cell_type"] == t))
      )
      # check if there are any unused cells of this type
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
        available.samples <- available.samples[-samples.to.add]
      }
    }

    # deconvolve with this training set
    temp.results <- deconvolute(
      training.expr = temp.expr,
      training.pheno = temp.pheno,
      test.expr = test.exprs,
      test.pheno = test.pheno,
      algorithms = algorithms,
      verbose = verbose,
      split.data = FALSE,
      exclude.from.signature = exclude.from.signature,
      optimize = TRUE,
      max.genes = NULL,
      n.bulks = 0,
      bulks = bulk.data,
      n.repeats = 1
    )
    sample.size.lists[[as.character(step*step.size)]][[as.character(rep)]] <- temp.results[[1]][[1]]
  }
}
sample.size.lists[["bulk.props"]] <- temp.results$bulk.props
return(sample.size.lists)
}

