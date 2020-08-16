library(DAB)

cll_data <- read_data("/data/tim/cll_data_uncompressed.h5")

cll.exprs <- cll_data$sc.counts
cll.pheno <- cll_data$sc.pheno

#cll.samples <- sample(1:ncol(cll.exprs), ceiling(0.25 * ncol(cll.exprs)), replace = F)
#cll.exprs <- cll.exprs[,cll.samples]
#cll.pheno <- cll.pheno[cll.samples,]

bulk.exprs <- cll_data$bulk.counts
bulk.props <- cll_data$bulk.props

# remove patient '12'; does not contain useful cells
to.remove <- which(nchar(as.character(cll.pheno$patient)) < 4)
cll.exprs <- cll.exprs[, -to.remove]
cll.pheno <- cll.pheno[-to.remove, ]

# split by patient
patient.info <- substr(as.character(cll.pheno$patient), 1, 1)

training.samples <- which(patient.info == 5 | patient.info == 1)
test.samples <- which(patient.info == 6 | patient.info == 8)
grouping <- rep(1, (length(training.samples)+length(test.samples)))
grouping[test.samples] <- 2

genesets <- readRDS("development/genesets.RDS")

cll.pheno <- cbind(cll.pheno, sample.name = colnames(cll.exprs))

benchmark(
  sc.counts = cll.exprs
  , sc.pheno = cll.pheno
  , real.counts = bulk.exprs
  , real.props = bulk.props
  , benchmark.name = "cll_benchmark"
  , exclude.from.signature = c("unassigned")
  , genesets = genesets
  , simulation.bulks = TRUE
  , simulation.genes = TRUE
  , simulation.samples = TRUE
  , simulation.subtypes = TRUE
  , repeats = 5
  , grouping = as.factor(grouping),
  verbose = TRUE,
  cell.type.column = "celltype",
  temp.dir = "/data/tim/DAB_complete/",
  input.algorithms = list("DTD", "Least_Squares", "MuSiC", "DeconRNASeq" ,"CIBERSORT", "BSEQ-sc")
)

tirosh_data <- read_data("/data/tim/tirosh_uncompressed_hgnc.h5")
tirosh.exprs <- tirosh_data$sc.counts
tirosh.pheno <- tirosh_data$sc.pheno

patients <- unique(tirosh.pheno$patient)
training.patients <- patients[sample(1:32, 20, replace = F)]
training.samples <- which(tirosh.pheno$patient %in% training.patients)

grouping <- rep(2, nrow(tirosh.pheno))
grouping[training.samples] <- 1

benchmark(
	  sc.counts = tirosh.exprs,
	  sc.pheno = tirosh.pheno,
	  real.counts = NULL,
	  real.props = NULL,
	  benchmark.name = "tirosh_benchmark",
	  exclude.from.signature = c("unknown"),
	  genesets = NULL,
	  simulation.bulks = TRUE,
	  simulation.genes = TRUE, 
	  simulation.samples = TRUE,
	  simulation.subtypes = TRUE,
	  repeats = 5,
	  grouping = as.factor(grouping),
	  verbose = TRUE,
	  temp.dir = "/data/tim/DAB_complete/",
	  sample.name.column = "sample_id",
	  input.algorithms = list("DTD", "Least_Squares", "MuSiC", "DeconRNASeq", "CIBERSORT", "BSEQ-sc")
)

# cll bulks with tirosh data

# match cell type labels
tirosh.pheno$cell_type[which(tirosh.pheno$cell_type == "T_CD4")] <- "cd4+"
tirosh.pheno$cell_type[which(tirosh.pheno$cell_type == "T_CD8")] <- "cd8+"
tirosh.pheno$cell_type[which(tirosh.pheno$cell_type == "macro")] <- "mono"
tirosh.pheno$cell_type[which(tirosh.pheno$cell_type == "NK")] <- "nk"

benchmark(
	  sc.counts = tirosh.exprs,
	  sc.pheno = tirosh.pheno,
	  real.counts = bulk.exprs,
	  real.props = bulk.props,
	  benchmark.name = "tirosh_on_cll",
	  genesets = NULL,
	  simulation.bulks = FALSE,
	  simulation.genes = FALSE,
	  simulation.samples = FALSE,
	  simulation.subtypes = FALSE,
	  repeats = 5,
	  grouping = factor(rep(0, ncol(tirosh.exprs)), levels = c(0,1)),
	  verbose = TRUE,
	  temp.dir = "/data/tim/DAB_complete",
	  sample.name.column = "sample_id",
	  input.algorithms = list("DTD", "Least_Squares", "DeconRNASeq", "MuSiC", "CIBERSORT", "BSEQ-sc")
)
