library(DeconvolutionAlgorithmBenchmarking)
# data from PB2
load("development/cll_normalized_sampled.rda")
# there are duplicated rows to deal with
cll.exprs <- cll.exprs[-which(duplicated(cll.exprs)), ]
bulks <- readRDS("development/real_bulks.rds")
genesets <- readRDS("development/genesets.RDS")
# remove patient '12'; does not contain useful cells
to.remove <- which(nchar(as.character(cll.pheno$patient)) < 4)
cll.exprs <- cll.exprs[, -to.remove]
 cll.pheno <- cll.pheno[-to.remove, ]

# split by patient
patient.info <- substr(as.character(cll.pheno$patient), 1, 1)
  
# patient.info
# 1    5    6    8
# 439 2035 1524  379
# try to make it even, but training set slightly larger:
# also, 5 and 8 are rather similar overall (dendrogram),
# so mix it for better performance
training.samples <- which(patient.info == 5 | patient.info == 1)
test.samples <- which(patient.info == 6 | patient.info == 8)
grouping <- rep(1, (length(training.samples)+length(test.samples)))
grouping[test.samples] <- 2

benchmark(cll.exprs, cll.pheno, bulks$bulks, bulks$props, benchmark.name = "test_benchmark", exclude.from.signature = c("unassigned", "not_annotated"), genesets = genesets, temp.dir = "/tmp/", simulation = c("genes" = T, "samples" = T, "bulks" = T), repeats = 3, grouping = as.factor(grouping))
