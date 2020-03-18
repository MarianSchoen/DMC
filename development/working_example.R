library(DeconvolutionAlgorithmBenchmarking)
# data from PB2
load("Git/WS1920-Mirus-DeconvAlgoComparison/data/cll_normalized_sampled.rda")
# there are duplicated rows to deal with
cll.exprs <- cll.exprs[-which(duplicated(cll.exprs)), ]
bulks <- readRDS("Git/WS1920-Mirus-DeconvAlgoComparison/data/real_bulks.rds")
genesets <- readRDS("Git/WS1920-Mirus-DeconvAlgoComparison/data/genesets.RDS")

benchmark(cll.exprs, cll.pheno, bulks$bulks, bulks$props, benchmark.name = "test_benchmark", exclude.from.signature = c("unassigned", "not_annotated"), genesets = genesets, temp.dir = "/data/temp/", simulation = c("genes" = T, "samples" = T, "bulks" = T), repeats = 3)
