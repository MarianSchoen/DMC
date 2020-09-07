#' main function of the deconvolution benchmark
#'
#' @param sc.counts non-negative numeric matrix with features as rows, and 
#' scRNA-Seq profiles as columns. \code{ncol(sc.counts)} must equal \code{nrow(sc.pheno)}
#' @param sc.pheno data frame with scRNA-Seq profiles as rows, and pheno entries
#'  in columns. \code{nrow(sc.pheno)} must equal \code{ncol(sc.counts)}. Cell types need 
#'  to be specified in a column named `cell.type.column` and the patient/origin
#'  (if available) in a column named 'patient'
#' @param bulk.counts non-negative numeric matrix, with features as rows, and 
#' bulk RNA-Seq profiles as columns. \code{ncol(sc.counts)} must equal 
#' \code{nrow(bulk.props)}
#' @param bulk.props non-negative numeric matrix specifying the amount of each
#' cell type in all each bulk, with cell types as rows and bulk RNA-Seq profiles as columns.
#' @param benchmark.name string, name of the benchmark. Will be used as name
#' for the results directory
#' @param grouping factor with 2 levels, and \code{length(grouping)} must be 
#' \code{ncol(sc.counts)}. Assigns each scRNA-Seq profile to either 
#' test or train cohort. 
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information?
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @param sample.name.column string, which column of 'pheno'
#' holds the sample name information; optional, default NULL
#' @param input.algorithms list containing a list for each algorithm. 
#' Each sublist contains \cr 1) name: character \cr 2) algorithm: function \cr
#' For predefined algorithms it is sufficient to supply only the name instead of the sublist,
#' e.g. \code{algorithms = list(list(name = 'DTD', algorithm = run_dtd), "DTD")}. \cr
#' If no list is supplied, all implemented algorithms 
#' (CIBERSORT, DeconRNASeq, DTD, Least_Squares, BSEQ-sc and MuSiC) are selected.
#' @param simulation.bulks boolean, should deconvolution of simulated bulks be
#' performed? default: FALSE
#' @param simulation.genes boolean, should deconvolution of simulated bulks with 
#' predefined genesets be performed? default: FALSE
#' @param simulation.samples boolean, should deconvolution of simulated bulks with
#' varying number of randomly selected training profiles be performed? default: FALSE
#' @param simulation.subtypes boolean, should deconvolution of simulated bulks with
#' artificial subtypes of given cell types be performed? default: FALSE
#' @param genesets list of string vector, must match 'rownames(sc.counts)'. default: NULL
#' @param metric method of result evaluation; either "cor" (default) or a function. An 
#' evaluation function must receive two vectors (real and estimated proportions of a cell type across all bulks)
#' as input and return a single value as output. Note that the score should be between 0 and 1
#' @param metric.name string, name of the evaluation metric used; not needed if metric is a string ("cor"). If metric is a function and metric.name
#' is NULL, the default will be "custom metric"
#' @param repeats numeric > 0, number of repetitions for each algorithm in each setting.
#' default: 3
#' @param temp.dir string, directory where data, and benchmarks get stored. default: NULL,
#' will be set to '.tmp'
#' @param exclude.from.bulks vector of strings, cell types that should not be 
#' predicted by the algorithms. default: NULL
#' @param exclude.from.signature vector of strings, cell types that should not be 
#' included in the simulated bulks. default: NULL
#' @param n.bulks numeric > 0, number of bulks to simulate. default 500
#' @param cpm boolean, should the sc profiles and bulks be scaled to counts
#' per million? default: TRUE
#' @param verbose boolean, should progress information be printed to the screen? default: FALSE
#' @param avg.profiles.per.subcluster integer vector specifying the average number of profiles assigned to each subtype in each simulation step. 
#' If not specified, a vector of length 'n.cluster.sizes'  will be automatically generated based on the data.
#' @param n.cluster.sizes integer, number of subtype cluster sizes to generate if avg.profiles.per.subcluster is not specified
#' @param cibersort.path string, path to CIBERSORT source code. Necesssary, if cibersort or bseq-sc wrappers are used.
#' @param n.profiles.per.bulk positive numeric, number of samples to be randomly
#' drawn for each simulated bulk; default 1000;
#' 
#' @return NULL, results are stored via hdf5 to 'temp.dir'
#' @export
#' 
#' @examples see either 'working_example.R', or 'working_example_fast.R'
benchmark <- function(
  sc.counts, 
  sc.pheno, 
  bulk.counts, 
  bulk.props,
  benchmark.name,
  grouping,
  cell.type.column = "cell_type",
  patient.column = "patient",
  sample.name.column = "sample.name",
  input.algorithms = NULL,
  simulation.bulks = FALSE,
  simulation.genes = FALSE,
  simulation.samples = FALSE,
  simulation.subtypes = FALSE,
  genesets = NULL, 
  metric = "cor",
  metric.name = NULL,
  repeats = 5, 
  temp.dir = NULL, 
  exclude.from.bulks = NULL, 
  exclude.from.signature = NULL, 
  n.bulks = 500, 
  cpm = TRUE, 
  verbose = FALSE,
  avg.profiles.per.subcluster = NULL,
  n.cluster.sizes = 5,
  cibersort.path = NULL,
  n.profiles.per.bulk = 1000
  ){
	  if(verbose) tictoc::tic("Benchmark")
	if(verbose) cat("calculating checksum\n")
	hash <- digest::digest(list(sc.counts, sc.pheno, bulk.counts, bulk.props, benchmark.name, grouping, exclude.from.bulks, exclude.from.signature, n.bulks, cpm, n.cluster.sizes, genesets, avg.profiles.per.subcluster))
	# check whether temporary directory is available and writeable
	# if not specified use .tmp in working directory
	if(is.null(temp.dir)){
		warning("No temporary directory was provided. Using .tmp in current directory.")
		temp.dir <- paste(getwd(),"/.tmp",sep = "")
	}
	# check that temp.dir exists and is writeable 
	# stop if it is not
	if(!dir.exists(temp.dir)){
		flag <- dir.create(temp.dir, recursive = TRUE)
		if(!flag) stop("Could not create temp directory. Please provide a writeable directory.")
	}else{
		if(is.null(benchmark.name) || benchmark.name == ""){
			stop("Invalid benchmark name. Please provide a unique name for this benchmark.")
		}else{
			if(dir.exists(paste(temp.dir, benchmark.name,sep="/"))){
				message("Found existing project directory within temp. Using present data where possible.")
				# compare to the old hash and proceed only if it does not exist or they are similar
				if(file.exists(paste(temp.dir, benchmark.name, "hash.rds",sep="/"))){
					hash.old <- readRDS(paste(temp.dir, benchmark.name, "hash.rds", sep="/"))
					if(hash != hash.old){
						stop("Hash values of current and old function call do not match. If you changed your data or other important parameters please name your benchmark differently.")
					}
				}
				saveRDS(hash, paste(temp.dir, benchmark.name, "hash.rds",sep="/"))
				rm("hash")
			}else{
				flag <- dir.create(paste(temp.dir, "/", benchmark.name, sep = ""))
				if(!flag) {
					stop("Could not create project folder. Is temp writeable?")
				}
			}
		}
	}
	output.folder <- paste(temp.dir,"/",benchmark.name,sep="")

	# check input parameters
	if(ncol(sc.counts) != nrow(sc.pheno)){
		stop("Dimensions of sc.counts and sc.pheno do not match")
	}
	if(!is.null(bulk.counts) && !is.null(bulk.props)){
	if(ncol(bulk.counts) != ncol(bulk.props)){
		stop("Number of bulks in bulk.counts and bulk.props do not match")
	}
	}
	if(is.null(genesets)){
		warning("No gene sets provided; skipping that benchmark")
		simulation.genes <- FALSE
	}
	if(is.character(metric)){
		if(metric != "cor"){
			stop("metric must be either \"cor\" or a function")
		}else{
			metric <- cor
			if(is.null(metric.name) || !is.character(metric.name)){
				metric.name <- "pearson correlation"
			}
		}
	}else{
		if(is.function(metric)){
			warning("Using custom evaluation function / metric. Unexpected results and plots may occur.")
			if(is.null(metric.name) || !is.character(metric.name)){
				metric.name <- "custom metric"
			}
		}else{
			stop("Function corresponding to 'metric' could not be found.")
		}
	}
	if(!is.factor(grouping) || !length(levels(grouping)) == 2 || !length(grouping) == ncol(sc.counts)){
		stop("Invalid sample grouping. Must be a factor of length nrow(sc.counts) with two levels indicating training and validation set")
	}
	if(!all(is.logical(c(simulation.bulks, simulation.genes, simulation.samples, simulation.subtypes)))){
		stop("Invalid value for at least one benchmark indicator. Have to be logical.")
	}
	if(!is.numeric(repeats)){
		stop("Invalid number of repeats. Must be numeric.")
	}
	if(!is.numeric(n.bulks)){
		stop("Invalid input. n.bulks must be numeric.")
	}
	if(!is.numeric(n.cluster.sizes)){
		stop("Invalid input. n.cluster.sizes must be integer")
	}
	if(repeats < 1){
		warning("Repeats must be greater than 0. Setting to 1.")
		repeats <- 1
	}
	if(n.cluster.sizes < 1 || as.integer(n.cluster.sizes) != n.cluster.sizes){
		stop("n.cluster.sizes must be integer > 0")
	}
	if(n.bulks < 1){
		warning("Number of bulks must be greater than 0. Setting to default.")
		n.bulks <- 1000
	}
	if(!is.null(exclude.from.bulks)){
		if(!all(exclude.from.bulks %in% unique(sc.pheno[[cell.type.column]]))){
			stop("Unknown cell type(s) in exclude.from.bulks. Please select only cell types present in pheno data.")
		}
	}
	if(!is.null(exclude.from.signature)){
		if(!all(exclude.from.signature %in% unique(sc.pheno[[cell.type.column]]))){
			stop("Unknown cell types(s) in exclude.from.signature. Please select only cell types present in pheno data.")
		}
	}
	if(!rmarkdown::pandoc_available()){
		warning("pandoc was not found on your system. you can perform the benchmark but will not be able to render the report until this dependency is fulfilled.")
	}
	if(!is.null(avg.profiles.per.subcluster)){
		if(!is.numeric(avg.profiles.per.subcluster)){
			stop("avg.profiles.per.subcluster must be a vector of integers")
		}
		if(any(as.integer(avg.profiles.per.subcluster) != avg.profiles.per.subcluster)){
			stop("avg.profiles.per.subcluster must be a vector of integers")
		}
	}

	use.cibersort <- FALSE
	if(exists("CIBERSORT") && is.function(CIBERSORT)){
		use.cibersort <- TRUE
	}else{
		if(!is.null(cibersort.path)){
			if(file.exists(cibersort.path)){
				source(cibersort.path)
				if(!exists("CIBERSORT")){
					warning("Could not source CIBERSORT.")
				}else{
					use.cibersort <- TRUE
				}
			}else{
				warning(paste("Specified file not found:", cibersort.path))
			}
		}else{
			warning("No CIBERSORT source file supplied. Any wrappers depending on CIBERSORT can not be used.")
		}
	}

	# check and process algorithms input
	algorithms <- list(
			   list(algorithm = run_dtd, name = "DTD"),
			   list(algorithm = run_cibersort, name = "CIBERSORT"),
			   list(algorithm = run_deconrnaseq, name = "DeconRNASeq"),
			   list(algorithm = run_least_squares, name = "Least_Squares"),
			   list(algorithm = run_music, name = "MuSiC"),
			   list(algorithm = run_bseqsc, name = "BSEQ-sc")
	)
	algorithm.names <- sapply(algorithms, function(x) x$name)

	if(!is.null(input.algorithms)){
		# input.algorithms must be a list
		if(!is.list(input.algorithms)){
			stop("Invalid algorithm input. input.algorithms must be a list.")
		}
		predef.algos <- c()
		new.algos <- list()
		# allow either lists containing two entries (for function and name)
		# or characters specifying one of the predefined algorithms
		for(a in input.algorithms) {
			if(is.list(a)){
				# for now check only whether algorithm exists, not its output
				if(exists(as.character(substitute(a$algorithm))) && is.character(a$name)){
					new.algos <- c(new.algos, list(a))
				}else{
					if(!is.character(a$name)){
						stop("Invalid algorithm.")
					}else{
						stop(paste0("Invalid algorithm: ", a$name))
					}
				}
			}else{
				if(is.character(a) && a %in% algorithm.names){
					predef.algos <- c(predef.algos, which(algorithm.names == a))	
				}else{
					stop("Invalid algorithm. User supplied algorithms can be given as a list with two entries: name (name of the algorithm) and algorithm (wrapper to call the algorithm)")
				}
			}
		}
		# join supplied (new) algorithms and specified existing algorithms
		if(length(predef.algos) > 0 || length(new.algos) > 0){
			if(length(predef.algos) > 0) {
				algorithms <- algorithms[predef.algos]
				if(length(new.algos) > 0) {
					algorithms <- c(algorithms, new.algos)
				}
			}else{
				algorithms <- new.algos
			}
		}else{
			stop("No algorithms selected")
		}
		algorithm.names <- sapply(algorithms, function(x) x$name)
	}

	if(("CIBERSORT" %in% algorithm.names || "BSEQ-sc" %in% algorithm.names) && !use.cibersort){
		warning("CIBERSORT and/or BSEQ-sc are among the specified algorithms but CIBERSORT is not available. Removing these algorithms...")
		to.use <- which(!algorithm.names %in% c("BSEQ-sc", "CIBERSORT"))
		algorithms <- algorithms[to.use]
		algorithm.names <- algorithm.names[to.use]
	}

	check_algorithms(algorithms)
	
	# Data preparation
	cat("data preparation...\t\t", as.character(Sys.time()), "\n", sep = "")

	# load / process / store data
	# if it exists load previously processed data from temp
	if(file.exists(paste(output.folder,"/input_data/training_set.h5",sep="")) && file.exists(paste(output.folder,"/input_data/validation_set.h5",sep=""))){
		if(verbose) cat("Loading data found in temp directory\n")
		training_set <- read_data(paste(output.folder, "/input_data/training_set.h5", sep = ""))
		validation_set <- read_data(paste(output.folder, "/input_data/validation_set.h5", sep=""))
		training.exprs <- training_set$sc.counts
		training.pheno <- training_set$sc.pheno
		test.exprs <- validation_set$sc.counts
		test.pheno <- validation_set$sc.pheno
		sim.bulks <- list(bulks = validation_set$bulk.counts, props = validation_set$bulk.props)
		
		# don't need these any more
		rm(list = c("training_set", "validation_set"))
	}else{
		dir.create(paste(output.folder,"/input_data",sep=""), recursive = TRUE)
	}
	# save input data of benchmark() to temp directory
	function.call <- match.call()
	if(!file.exists(paste(output.folder, "input_data/raw.h5", sep = "/"))){
	  suppressMessages(suppressWarnings(write_data(sc.counts, sc.pheno, bulk.counts, bulk.props, filename = paste(output.folder,"input_data/raw.h5", sep = "/"))))
	}
	if(!file.exists(paste(output.folder, "input_data/params.h5", sep = "/"))){
	  suppressMessages(suppressWarnings(write_misc_input(algorithm.names = algorithm.names, genesets = genesets, function.call = function.call, grouping = grouping, file = paste(output.folder,"input_data/params.h5",sep="/"))))
	}else{
	  # join list of algorithms and genesets from previous function calls
	  input_params <- read_misc_input(paste(output.folder, "input_data/params.h5", sep = "/"))
	  suppressMessages(suppressWarnings(write_misc_input(algorithm.names = unique(c(algorithm.names, input_params$algorithms)), genesets = unique(c(genesets, input_params$genesets)), function.call = function.call, grouping = grouping, file = paste(output.folder,"input_data/params.h5",sep="/"))))
	  # don't need this any more
	  rm("input_params")
	}
	
	# free ram
	gc()

	# if any of the required data is missing preprocess input data for deconvolution
	if(!exists("training.exprs") || !exists("training.pheno") || !exists("test.exprs") || !exists("test.pheno") || !exists("sim.bulks")){
		if(verbose) cat("preprocessing data\n")
		if(cpm){
			if(verbose) cat("scaling expression profiles to cpm\n")
			sc.counts <- scale_to_count(sc.counts)
			gc()
			if(!is.null(bulk.counts))
				bulk.counts <- scale_to_count(bulk.counts)
		}

		if(verbose) cat("splitting into test and training data\n")
		# split data into test and validation set
		if(length(unique(grouping)) == 2){
			split.data <- split_dataset(sc.counts, sc.pheno, method = "predefined", grouping = grouping)
			training.exprs <- split.data$training$exprs
			training.pheno <- split.data$training$pheno
			test.exprs <- split.data$test$exprs
			test.pheno <- split.data$test$pheno
			
			# remove unnecessary variable
			rm("split.data")
		}else if(length(unique(grouping)) == 1){
			message("Grouping vector contains only one group. Disabling all simulations.")
			training.exprs <- sc.counts
			training.pheno <- sc.pheno
			test.exprs <- NULL
			test.pheno <- NULL
			simulation.bulks = FALSE
			simulation.genes = FALSE
			simulation.samples = FALSE
			simulation.subtypes = FALSE
		}else{
			stop("Invalid grouping vector")
		}
	  
	  # this is not needed any more
	  rm(list = c("sc.counts", "sc.pheno"))
	  gc()

		# exclude.from.bulks is a parameter of benchmark(), create_bulks expects the opposite
		if(!is.null(exclude.from.bulks) && length(intersect(unique(test.pheno[[cell.type.column]]), exclude.from.bulks)) > 0){
			include.in.bulks <- unique(test.pheno[[cell.type.column]])[-which(unique(test.pheno[[cell.type.column]] %in% exclude.from.bulks))]
		}else{
			include.in.bulks <- NULL
		}

		# create simulated bulks if test data is available and they are needed
	  # they are only needed if bulk/geneset/training set simulations are performed
	  # subtype simulation generates its own bulks
		if(!is.null(test.exprs) && any(c(simulation.bulks, simulation.genes, simulation.samples))){
			cat("creating artificial bulks for simulation\n")
			sim.bulks <- create_bulks(
			  test.exprs, 
			  test.pheno, 
			  n.bulks, 
			  include.in.bulks, 
			  sum.to.count = cpm, 
			  cell.type.column = cell.type.column, 
			  n.profiles.per.bulk = n.profiles.per.bulk
			)
		}else{
			sim.bulks <- list(bulks = NULL, props = NULL)
		}
		# save processed data for later use and repeated benchmarks
		suppressMessages(suppressWarnings(write_data(training.exprs, training.pheno, filename = paste(output.folder, "/input_data/training_set.h5", sep = ""))))
		suppressMessages(suppressWarnings(write_data(test.exprs, test.pheno, sim.bulks$bulks, sim.bulks$props, paste(output.folder, "/input_data/validation_set.h5", sep=""))))
	}

	# rescale real quantities to be between 0 and 1
	if(!is.null(bulk.props)){
		bulk.props <- apply(bulk.props, 2, function(x) {x / sum(x)})
	}

	# assume that samples in expression and pheno data are in the correct order
	# if names do not match, assign sample names from expression to pheno data
	if(!identical(rownames(training.pheno), colnames(training.exprs))){
		rownames(training.pheno) <- colnames(training.exprs)
		colnames(training.exprs) <- rownames(training.pheno)
	}
	if(!is.null(test.exprs)){
	  if(!identical(rownames(test.pheno), colnames(test.exprs))){
	    rownames(test.pheno) <- colnames(test.exprs)
	    colnames(test.exprs) <- rownames(test.pheno)
	  }
	}

	# make sure subtypes of types in exclude.from.signature are also excluded
	if(any(training.pheno[[cell.type.column]] %in% exclude.from.signature) && "subtype" %in% colnames(training.pheno)){
	  exclude.from.signature <- c(exclude.from.signature, unique(paste(training.pheno[[cell.type.column]] ,training.pheno$subtype, sep = ".")[which(training.pheno[[cell.type.column]] %in% exclude.from.signature)]))
	}
	if(any(test.pheno[[cell.type.column]] %in% exclude.from.signature) && "subtype" %in% colnames(test.pheno)){
	  exclude.from.signature <- c(exclude.from.signature, unique(paste(test.pheno[[cell.type.column]] ,test.pheno$subtype, sep = ".")[which(test.pheno[[cell.type.column]] %in% exclude.from.signature)]))
	}
	exclude.from.signature <- unique(exclude.from.signature)
	
	# Deconvolution
	cat("deconvolution...\t\t", as.character(Sys.time()), "\n", sep = "")

	# select algorithms that have not been evaluated in a previous run
	previous.results <- list()
	# read available files for this benchmark
	if(dir.exists(paste(output.folder,"/results/real/",sep=""))){
		files <- list.files(paste(output.folder, "/results/real/", sep = ""), full.names = T, pattern = "deconv.*.h5")
		if(length(files) > 0){
			for(i in 1:length(files)){
				f <- files[i]
				previous.results[[i]] <- read_result_list(f)
			}
		}
	}else{
		dir.create(paste(output.folder,"/results/real/",sep=""), recursive = T)
	}

	# exclude algorithms for which results are already present
	if(length(previous.results) == 0){
		to.run <- 1:length(algorithms)
	}else{
		present.algorithms <- c()
		for(r in previous.results){
			present.algorithms <- c(present.algorithms, unique(as.character(prepare_data(r)$algorithm)))
		}
		if(verbose) print(present.algorithms)
		to.run <- which(! algorithm.names %in% present.algorithms)
		rm(list = c("previous.results", "present.algorithms"))
		gc()
	}
	
	# remove variables
	
	
	res.no <- length(previous.results) + 1
	if(!is.null(bulk.counts) && !is.null(bulk.props)){
	# deconvolute real bulks
	cat("deconvolve real bulks...\t", as.character(Sys.time()), "\n", sep = "")
	if(length(to.run)>0){
		if(verbose) cat("algorithms:", sapply(algorithms[to.run], function(x) {x$name}), "\n", sep = " ")
		real.benchmark <- deconvolute(
			training.exprs, 
			training.pheno, 
			NULL, NULL, # do not need test data to deconvolve real bulks
			algorithms[to.run], 
			verbose, 
			TRUE, 
			NULL, 
			exclude.from.signature, 
			TRUE, 
			NULL, 
			0, 
			list(bulks = bulk.counts, props = bulk.props), 
			repeats,
			cell.type.column = cell.type.column,
			patient.column = patient.column
		)
		suppressMessages(suppressWarnings(write_result_list(real.benchmark, paste(output.folder, "/results/real/deconv_output_",res.no,".h5",sep=""))))

		# perform bootstrapping
		cat("bootstrapping...\t\t", as.character(Sys.time()), "\n", sep = "")
		estimates <- list()
		# use the results of the algorithms in the first repetition
		for(a in algorithms[to.run]){
			estimates[[a$name]] <- real.benchmark$results.list[["1"]][[a$name]]$est.props
		}
		props <- list(real = bulk.props, est = estimates)
		bootstrap.real <- bootstrap_bulks(props)
		h5_write_mat(bootstrap.real, paste(output.folder, "/results/real/bootstrap_bulks",res.no,".h5",sep=""))
	}
	}
	
	# remove real bulk results from RAM
	rm(list = c("real.benchmark", "estimates", "props", "bootstrap.real"))
	gc()

	# iterate through supplied simulation vector and perform those that are TRUE
	available.sims <- c(simulation.genes, simulation.samples, simulation.bulks, simulation.subtypes)
	names(available.sims) <- c("genes", "samples", "bulks", "subtypes")
	if(any(available.sims)) cat("starting simulations...\t\t", as.character(Sys.time()), "\n", sep = "")
	for(s in names(available.sims)){
		if(!available.sims[s]) next

		# read previous results and exclude present algorithms
		previous.results <- list()
		if(dir.exists(paste(output.folder, "/results/simulation/",s,"/", sep=""))){
			files <- list.files(paste(output.folder, "/results/simulation/",s,"/",sep = ""), full.names = T, pattern = "*.h5")
			if(length(files) > 0){
				for(i in 1:length(files)) {
					f <- files[i]
					previous.results[[i]] <- read_result_list(f)
				}
			}
		}else{
			dir.create(paste(output.folder,"/results/simulation/",s,"/",sep=""), recursive = T)
		}
		if(length(previous.results) == 0){
			to.run <- 1:length(algorithms)
		}else{
			present.algorithms <- c()
			for(r in previous.results){
				present.algorithms <- c(present.algorithms, unique(as.character(prepare_data(r)$algorithm)))
			}
			if(verbose) print(present.algorithms)
			to.run <- which(! algorithm.names %in% present.algorithms)
		}
		res.no <- length(previous.results) + 1
		
		# remove variables containing previous results
		rm(list = c("previous.results", "present.algorithms"))
		gc()

		# execute benchmark corresponding to s and save results
		if(length(to.run)>0){
			if(s == "bulks"){
				cat("bulk simulation...\t\t", as.character(Sys.time()), "\n", sep = "")
				benchmark.results <- deconvolute(training.exprs, training.pheno, NULL, NULL, algorithms[to.run], verbose, FALSE, NULL, exclude.from.signature, TRUE, NULL, 0, sim.bulks, repeats, cell.type.column = cell.type.column, patient.column = patient.column)
			}
			if(s == "genes"){
				cat("geneset simulation...\t\t", as.character(Sys.time()), "\n", sep = "")
			  benchmark.results <- geneset_benchmark(training.exprs, training.pheno, NULL, NULL, genesets, algorithms[to.run], sim.bulks, repeats, exclude.from.signature, verbose, split.data = FALSE, cell.type.column = cell.type.column, patient.column = patient.column)
			}
			if(s == "samples"){
				cat("sample simulation...\t\t", as.character(Sys.time()), "\n", sep = "")
			  benchmark.results <- sample_size_benchmark(training.exprs, training.pheno, NULL, NULL, algorithms[to.run], sim.bulks, repeats, exclude.from.signature, 0.25, verbose, split.data = FALSE, cell.type.column = cell.type.column, patient.column = patient.column)
			}
			if(s == "subtypes"){
				cat("subtype simulation...\t\t", as.character(Sys.time()), "\n", sep = "")
				if(!is.null(test.exprs) && !is.null(test.pheno)){
					all.exprs <- cbind(training.exprs, test.exprs)
					all.pheno <- rbind(training.pheno, test.pheno)
				}else{
					all.exprs <- training.exprs
					all.pheno <- training.pheno
				}
				
			  benchmark.results <- fine_coarse_subtype_benchmark(
					all.exprs, all.pheno, 
					cell.type.column = cell.type.column, 
					sample.name.column = sample.name.column, 
					verbose = verbose, 
					algorithm.list = algorithms[to.run],
					avg.profiles.per.subcluster = avg.profiles.per.subcluster, 
					n.cluster.sizes = n.cluster.sizes,
					patient.column = patient.column,
					n.bulks = n.bulks
				)
			  
			  rm(list = c("all.exprs", "all.pheno"))
			}
		  
			if(!dir.exists(paste(output.folder, "/results/simulation/", s, sep = ""))){
		 		dir.create(paste(output.folder, "/results/simulation/", s, sep = ""), recursive = TRUE)
			}
		  suppressMessages(suppressWarnings(write_result_list(benchmark.results, paste(output.folder, "/results/simulation/", s, "/deconv_output_", res.no, ".h5", sep = ""))))
		  rm("benchmark.results")
		}
	}

	cat("preparing results...\t\t", as.character(Sys.time()), "\n", sep = "")
	# deconvolution step is over
	# create results markdown
	
	report.path <- suppressWarnings(render_results(output.folder, metric, metric.name, benchmark.name, cell.type.column))
	cat("Done\t\t\t\t", as.character(Sys.time()), "\n\n", sep = "")
	cat("Report generated: ", report.path, "\n", sep = "")
	cat("Created plots can be found in: ", output.folder, "/report_plots/\n", sep = "")
	
	if(dir.exists("CIBERSORT")){
	  unlink("CIBERSORT", recursive = TRUE)
	}
	
	if(verbose) {
		t <- tictoc::toc(quiet = TRUE)
		time <- t$toc - t$tic
		if(time > 3600) {
			time <- time / 3600
			cat(t$msg, ": ", as.integer(time), "h ", as.integer((time-as.integer(time))*60), "m elapsed\n", sep = "")
		}else{
			if(time > 60) {
				time <- time / 60
				cat(t$msg, ": ", as.integer(time), "m ", as.integer((time-as.integer(time))*60), "s elapsed\n", sep = "")
			} else {
				cat(t$msg, ": ", time, "s elapsed\n", sep = "")
			}
		}
	}
}
