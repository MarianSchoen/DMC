#' main function of the benchmark
#'
#' @param sc.counts numeric matrix with features as rows, and scRNA-Seq profiles
#' as columns. 
#' @param sc.pheno 
#' @param real.counts 
#' @param real.props 
#' @param benchmark.name
#' @param grouping
#' @param input.algorithms 
#' @param simulations 
#' @param genesets 
#' @param metric 
#' @param repeats 
#' @param temp.dir 
#' @param exclude.from.bulks 
#' @param exclude.from.signature
#' @param n.bulks
#' @param cpm
#'
#' @return
#' @export
#'
#' @examples
benchmark <- function(sc.counts, 
					  sc.pheno, 
					  real.counts, 
					  real.props,  
					  benchmark.name,
					  grouping,
					  input.algorithms = NULL, 
					  simulations=c("bulks"=TRUE, "genes"=TRUE, "samples"=TRUE), 
					  genesets = NULL, 
					  metric = "cor", 
					  repeats = 3, 
					  temp.dir = NULL, 
					  exclude.from.bulks = NULL, 
					  exclude.from.signature = NULL, 
					  n.bulks = 1000, 
					  cpm = TRUE, 
					  verbose = FALSE){
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
			if(dir.exists(paste(temp.dir,"/",benchmark.name,sep=""))){
				print("Found existing project directory within temp")
			}else{
				flag <- dir.create(paste(temp.dir, "/", benchmark.name, sep = ""))
				if(!flag) {
					stop("Could not create project folder. Is temp writeable?")
				}
			}
		}
	}
	
	output.folder <- paste(temp.dir,"/",benchmark.name,sep="")

	# check counts and props input
	if(ncol(sc.counts) != nrow(sc.pheno)){
		stop("Dimensions of sc.counts and sc.pheno do not match")
	}
	if(ncol(real.counts) != ncol(real.props)){
		stop("Number of bulks in real.counts and real.props do not match")
	}
	if(is.null(genesets)){
		warning("No gene sets provided skipping this benchmark")
		simulations["genes"] <- F
	}
	if(!metric %in% c("cor", "mad", "rmsd")){
		stop("metric must be one of 'cor', 'mad', 'rmsd'")
	}
	if(!is.factor(grouping) || !length(levels(grouping)) == 2 || !length(grouping) == ncol(sc.counts)){
		stop("Invalid sample grouping. Must be a factor of length nrow(sc.counts) with two levels indicating training and validation set")
	}

	# check and process algorithms input
	algorithms <- list(
			   list(algorithm = run_dtd, name = "DTD"),
			   list(algorithm = run_cibersort, name = "CIBERSORT"),
			   list(algorithm = run_deconrnaseq, name = "DeconRNASeq"),
			   list(algorithm = run_least_squares, name = "Least_Squares")
	)
	algorithm.names <- sapply(algorithms, function(x) x$name)

	if(!is.null(input.algorithms)){
		# input.algorithms must be a list
		if(!is.list(input.algorithms)){
			stop("Invalid algorithm input")
		}
		predef.algos <- c()
		new.algos <- list()
		# allow either lists containing two entries (for function and name)
		# or characters specifying one of the predefined algorithms
		for(a in input.algorithms) {
			if(is.list(a)){
				# for now check only whether algorithm exists, not its output
				if(exists(as.character(substitute(a$algorithm))) && is.character(a$name)){
					new.algos <- c(new.algos, a)
				}else{
					stop("Invalid algorithm")
				}
			}else{
				if(is.character(a) && a %in% algorithm.names){
					predef.algos <- c(predef.algos, which(algorithm.names == a))	
				}else{
					stop("Invalid algorithm")
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


	# load / process / store data
	# if it exists load previously processed data from temp
	if(file.exists(paste(output.folder,"/input_data/processed.rda",sep=""))){
		if(verbose) print("Using data found in temp directory")
		load(paste(output.folder,"/input_data/processed.rda",sep=""))
	}else{
		dir.create(paste(output.folder,"/input_data",sep=""))
	}
	# save input data of benchmark() to temp directory
	function.call <- match.call()
	save(sc.counts, sc.pheno, real.counts, real.props, algorithm.names, genesets, function.call, file = paste(output.folder,"/input_data/raw.rda",sep=""))

	# if any of the required data is missing preprocess input data for deconvolution
	if(!exists("training.exprs") || !exists("training.pheno") || !exists("test.exprs") || !exists("test.pheno") || !exists("sim.bulks")){
		if(cpm){
			sc.counts <- scale_to_count(sc.counts)
			real.counts <- scale_to_count(real.counts)
		}
		# split (randomly at the moment) into test and validation set
		split.data <- split_dataset(sc.counts, sc.pheno, method = "predefined", grouping = grouping)
		training.exprs <- split.data$training$exprs
		training.pheno <- split.data$training$pheno
		test.exprs <- split.data$test$exprs
		test.pheno <- split.data$test$pheno
		# exclude.from.bulks is a paramweter of benchmark(), create_bulks expects the opposite
		if(!is.null(exclude.from.bulks) && length(intersect(unique(test.pheno$cell_type), exclude.from.bulks)) > 0){
			include.in.bulks <- unique(test.pheno$cell_type)[-which(unique(test.pheno$cell_type %in% exclude.from.bulks))]
		}else{
			include.in.bulks <- NULL
		}
		sim.bulks <- create_bulks(test.exprs, test.pheno, n.bulks, include.in.bulks, sum.to.count = cpm)
		# save processed data for later use and repeated benchmarks
		save(training.exprs, training.pheno, test.exprs, test.pheno, sim.bulks, file = paste(output.folder, "/input_data/processed.rda", sep=""))	
	}
	
	# we have not agreed on whether data plots should be generated yet ...
	#
	#
	#
	# begin deconvolution part

	# select algorithms that have not been evaluated in a previous runs
	if(verbose) print("Benchmarking performance on real bulks")
	previous.results <- list()
	# read available files for this benchmark
	if(dir.exists(paste(output.folder,"/results/real/",sep=""))){
		files <- list.files(paste(output.folder, "/results/real/", sep = ""), full.names = T, pattern = "*.rds")
		if(length(files) > 0){
			for(i in 1:length(files)){
				f <- files[i]
				previous.results[[i]] <- readRDS(f)
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
	}
	res.no <- length(previous.results) + 1

	# deconvolute real bulks
	if(length(to.run)>0){
		real.benchmark <- deconvolute(training.exprs, training.pheno, NULL, NULL, algorithms[to.run], verbose, FALSE, NULL, exclude.from.signature, TRUE, NULL, 0, list(bulks = real.counts, props = real.props), repeats)
		saveRDS(real.benchmark, paste(output.folder, "/results/real/deconv_output_",res.no,".rds",sep=""))
	}

	# iterate through supplied simulation vector and perform those that are TRUE
	available.sims <- c("genes", "samples", "bulks")
	if(any(simulations) && verbose) print("Starting simulations")
	for(s in names(simulations)){
		if(!simulations[s] || ! s %in% available.sims) next
		# read previous results and exclude present algorithms
		previous.results <- list()
		if(dir.exists(paste(output.folder, "/results/simulation/",s,"/", sep=""))){
			files <- list.files(paste(output.folder, "/results/simulation/",s,"/",sep = ""), full.names = T, pattern = "*.rds")
			if(length(files) > 0){
			for(i in 1:length(files)) {
				f <- files[i]
				previous.results[[i]] <- readRDS(f)
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

		# create unified interface for the benchmarks in the future?
		# execute benchmark corresponding to s and save results
		if(length(to.run)>0){
			if(s == "bulks"){
				print("bulk simulation")
				sim.bulk.benchmark <- deconvolute(training.exprs, training.pheno, NULL, NULL, algorithms[to.run], verbose, FALSE, NULL, exclude.from.signature, TRUE, NULL, 0, sim.bulks, repeats)
				benchmark.results <- sim.bulk.benchmark
			}
			if(s == "genes"){
				print("geneset simulation")
				sim.genes.benchmark <- geneset_benchmark(training.exprs, training.pheno, NULL, NULL, genesets, algorithms[to.run], sim.bulks, repeats, exclude.from.signature, verbose)
				benchmark.results <- sim.genes.benchmark
			}
			if(s == "samples"){
				print("sample simulation")
				sim.sample.benchmark <- sample_size_benchmark(training.exprs, training.pheno, NULL, NULL, algorithms[to.run], bulks, repeats, exclude.from.signature, 0.25, verbose)
				benchmark.results <- sim.sample.benchmark
			}
			if(!dir.exists(paste(output.folder, "/results/simulation/", s, sep = ""))){
		 		dir.create(paste(output.folder, "/results/simulation/", s, sep = ""), recursive = TRUE)
			}
			saveRDS(benchmark.results, paste(output.folder, "/results/simulation/",s,"/deconv_output_",res.no,".rds",sep=""))
		}
	}

	# deconvolution step is over
	# create results markdown
	render_results(output.folder, metric)
}
