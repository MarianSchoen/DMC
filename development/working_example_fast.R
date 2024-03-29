# no filthy workspace: 
rm(list = ls())

library(tictoc)

# while developing, maybe devtools::load_all is better than library, 
# as current changes get loaded, and must not be installed
# library(DeconvolutionAlgorithmBenchmarking)
devtools::load_all(".")

source("development/via_dtd_get_random_data.R")
data.list <- via_dtd_get_random_data()

genesets <- list(
  "1" = sample(rownames(data.list$sc.data), size = 0.5*nrow(data.list$sc.data)),
  "2" = sample(rownames(data.list$sc.data), size = 0.5*nrow(data.list$sc.data))
)

grouping <- sample(
  x = c(1,2)
  , size = ncol(data.list$sc.data)
  , replace = TRUE
)
tmp.dir <- paste(getwd(),".tmp",sep="/")
unlink(
  x = tmp.dir
  , recursive = TRUE
  )


tic("fast benchmark")
benchmark(
  sc.counts = data.list$sc.data
  , sc.pheno = data.frame(
    "cell_type" = data.list$sc.pheno, 
    "names" = names(data.list$sc.pheno))
  , real.counts = data.list$bulks
  , real.props = data.list$bulks.pheno
  , benchmark.name = "fast_test_benchmark"
  , exclude.from.signature = c()
  , genesets = genesets
  , simulation.bulks = TRUE
  , simulation.samples = TRUE
  , simulation.subtypes = TRUE
  , repeats = 3
  , grouping = as.factor(grouping)
  , temp.dir = tmp.dir
  , input.algorithms = list("DTD", "Least_Squares")
  , verbose = FALSE
)
toc()

