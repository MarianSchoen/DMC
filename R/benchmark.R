#' main function of the deconvolution benchmark
#'
#' @param sc.counts non-negative numeric matrix with features as rows, and
#' scRNA-Seq profiles as columns. \code{ncol(sc.counts)} must equal
#' \code{nrow(sc.pheno)}. May also be sparse matrix (class 'dgCMatrix')
#' @param sc.pheno data frame with scRNA-Seq profiles as rows, and pheno entries
#'  in columns. \code{nrow(sc.pheno)} must equal \code{ncol(sc.counts)}.
#'  Cell types need to be specified in column `cell.type.column`,
#'  the patient/origin (if available) in column `patient.column` and
#'  the sample names in column `sample.name.column`
#' @param bulk.counts non-negative numeric matrix, with features as rows, and
#' bulk RNA-Seq profiles as columns. \code{ncol(sc.counts)} must equal
#' \code{nrow(bulk.props)}. May also be sparse matrix (class 'dgCMatrix')
#' @param bulk.props non-negative numeric matrix specifying the amount of each
#' cell type in all each bulk, with cell types as rows and
#' bulk RNA-Seq profiles as columns.
#' @param benchmark.name string, name of the benchmark. Will be used as name
#' for the results directory
#' @param grouping Can be either\cr
#' 1) factor with 2 levels, and \code{length(grouping)} must be
#' \code{ncol(sc.counts)}. Assigns each scRNA-Seq profile to either
#' test or train cohort. 1 marks training samples, 2 marks test samples.\cr
#' 2) string, name of a column in sc.pheno containing "training" or "test" for all cells
#' @param cell.type.column string, which column of 'sc.pheno'
#' holds the cell type information? default 'cell_type'
#' @param patient.column string, which column of 'sc.pheno'
#' holds the patient information; optional, default 'patient'
#' @param sample.name.column string, which column of 'sc.pheno'
#' holds the sample name information; optional, default 'sample.name'
#' @param input.algorithms list containing a list for each algorithm.
#' Each sublist contains \cr 1) name: character \cr 2) algorithm: function \cr
#' 3) model: model to be supplied to the algorithm, optional \cr
#' For predefined algorithms it is sufficient to supply only the name instead of
#' the sublist, e.g. \code{algorithms = list(list(name = 'DTD',
#' algorithm = run_dtd), "MuSiC")}. \cr
#' If no list is supplied (default), all implemented algorithms
#' (CIBERSORT, DeconRNASeq, DTD, Least_Squares, BSEQ-sc and MuSiC) are selected.
#' @param simulation.bulks boolean, should deconvolution of simulated bulks be
#' performed? default: FALSE
#' @param simulation.genes boolean, should deconvolution of simulated bulks with
#' predefined genesets be performed? default: FALSE
#' @param simulation.samples boolean,
#' should deconvolution of simulated bulks with varying number of randomly
#' selected training profiles be performed? default: FALSE
#' @param simulation.subtypes boolean, should deconvolution of simulated bulks
#' with artificial subtypes of given cell types be performed? default: FALSE
#' @param genesets named list of string vectors, each must match subset of
#' 'rownames(sc.counts)'. default: NULL
#' @param repeats numeric > 0,
#' number of repetitions for each algorithm in each setting. default: 5
#' @param temp.dir string, directory where data, and benchmarks get stored.
#' default: NULL, using directory '.tmp' in working directory
#' @param exclude.from.signature vector of strings, cell types that should not
#' be predicted by the algorithms. default: NULL
#' @param exclude.from.bulks vector of strings,
#' cell types that should not be included in the simulated bulks. default: NULL
#' @param n.bulks numeric > 0, number of bulks to simulate. default 500
#' @param cpm boolean, should the sc profiles and bulks be scaled to counts
#' per million? default: FALSE
#' @param verbose boolean, should progress information be printed to the screen?
#'  default: FALSE
#' @param n.cluster.sizes vector of integers, number of artificial subtypes
#'  to generate per cell type; default: c(1, 2, 4, 8)
#' @param n.profiles.per.bulk positive numeric, number of samples to be randomly,
#'  default: 1000
#' @param report boolean, should an HTML report be generated? deafult TRUE
#'
#' @return list of\cr
#' 1) report_path: report path (string), NULL if no report is generated\cr
#' 2) bulk_results: deconvolution results for real bulks, NULL if no real bulks were supplied
#' @export

benchmark <- function(
  sc.counts, sc.pheno, 
  bulk.counts, bulk.props, 
  benchmark.name, grouping,
  cell.type.column = "cell_type",
  patient.column = "patient",
  sample.name.column = "sample.name",
  input.algorithms = NULL,
  simulation.bulks = FALSE, 
  simulation.genes = FALSE,
  simulation.samples = FALSE, 
  simulation.subtypes = FALSE,
  genesets = NULL,
  repeats = 5,
  temp.dir = NULL,
  exclude.from.bulks = NULL,
  exclude.from.signature = NULL,
  n.bulks = 500,
  cpm = FALSE,
  verbose = FALSE,
  n.cluster.sizes = c(1, 2, 4, 8),
  n.profiles.per.bulk = 1000,
  report = TRUE
  ) {
  # load Matrix library once again to ensure proper data handling
  suppressWarnings(suppressMessages(library(Matrix, quietly = TRUE)))

  # measure time
  if (verbose) tictoc::tic("Benchmark")
	if (verbose) cat("calculating checksum\n")
  hash <- digest::digest(
    list(
      sc.counts, sc.pheno, bulk.counts, bulk.props,
      benchmark.name, grouping, exclude.from.bulks,
      exclude.from.signature, n.bulks, cpm,
      n.cluster.sizes, genesets
    )
  )

	# if temp.dir not specified, use .tmp in working directory
	if (is.null(temp.dir)) {
		warning(
  		"No temporary directory was provided.
  		 Using .tmp in current directory."
		)
		temp.dir <- paste(getwd(), "/.tmp", sep = "")
	}

	# check that temp.dir exists and is writeable
	if(!dir.exists(temp.dir)){
		flag <- dir.create(temp.dir, recursive = TRUE)
		if (!flag) {
		  stop(
        "Could not create temp directory. Please provide a writeable directory."
      )
		}
		flag <- dir.create(paste(temp.dir, benchmark.name, sep = "/"), recursive = TRUE)
    if (!flag) {
      stop(
        "Could not create benchmark directory.
        Please provide writeable location."
      )
    }
	}else{
		if (is.null(benchmark.name) || benchmark.name == "") {
			stop("Invalid benchmark name.
			     Please provide a unique name for this benchmark.")
		}else{
		  # if directory already exists, check hash and try to load data
			if (dir.exists(paste(temp.dir, benchmark.name, sep = "/"))) {
				message("Found existing project directory within temp.
				        Using present data where possible.")
				# compare to the old hash and proceed only if they are identical
				if (file.exists(paste(temp.dir, benchmark.name, "hash.rds", sep = "/"))) {
					hash.old <- readRDS(paste(
					  temp.dir, benchmark.name, "hash.rds", sep = "/"
					))
					if(hash != hash.old){
						stop("Hash values of current and old function call do not match.
						     If you changed your data or other important parameters
						     name your benchmark differently.")
					}
				}
			}else{
				flag <- dir.create(paste(temp.dir, "/", benchmark.name, sep = ""))
				if (!flag) {
					stop("Could not create project folder. Is temp writeable?")
				}
			}
		}
	}
	saveRDS(hash, paste(temp.dir, benchmark.name, "hash.rds",sep = "/"))
	rm("hash")
	output.folder <- paste(temp.dir, "/", benchmark.name, sep = "")

	# check input parameters
	# single-cell counts and pheno data need to match
	if (!is.numeric(bulk.props))
	{
		stop("bulk.props must be a numeric matrix")
	}
	if (length(class(bulk.props)) != 2 || class(bulk.props)[1] != "matrix")
	{
		stop("bulk.props must be a numeric matrix")
	}
	if (!is.null(sc.counts) || !is.null(sc.pheno)) {
  	if (ncol(sc.counts) != nrow(sc.pheno)) {
  		stop("Dimensions of sc.counts and sc.pheno do not match")
  	}
  	# remove empty profiles
  	if (any(Matrix::colSums(sc.counts) == 0)) {
  		to.remove <- which(Matrix::colSums(sc.counts) == 0)
  		sc.pheno <- sc.pheno[-to.remove, ]
  		sc.counts <- sc.counts[, -to.remove]
  	}
	} else {
	  grouping <- NULL
	}
	
	# bulk counts and pheno data need to match and cell types must
	# overlap in single cell data and bulks if a report is to be generated
	if (!is.null(bulk.counts) && !is.null(bulk.props)) {
  	if (ncol(bulk.counts) != ncol(bulk.props)) {
  		stop("Number of bulks in bulk.counts and bulk.props do not match")
  	}
  	if (length(intersect(
          rownames(bulk.props), unique(sc.pheno[[cell.type.column]])
        )) == 0 && report && !is.null(sc.pheno)
      ) {
  		stop("RNA-Seq bulk data was supplied but there is no overlap between
        cell types in single-cell data and bulk.props.
        If you want to run the deconvolution anyway,
        set 'report = FALSE'.
        Bulks will be deconvoluted, but no automated report will be generated."
      )
  	}
	}
	if (is.null(genesets)) {
		warning("No gene sets provided; skipping that benchmark")
		simulation.genes <- FALSE
	}
	
	if (!is.null(grouping) && (!is.null(sc.counts) || !is.null(sc.pheno))) {
  	if (!is.factor(grouping) || length(levels(grouping)) != 2 ||
  	   length(grouping) != ncol(sc.counts)) {
		if (is.character(grouping) && length(grouping) == 1 && grouping %in% colnames(sc.pheno))
		{
			if (length(unique(sc.pheno[[grouping]])) == 2 && all(c("training", "test") %in% sc.pheno[[grouping]]))
			{
				grouping <- factor(ifelse(sc.pheno[[grouping]] == "training", 1, 2), levels = c(1,2))
			} else {
				stop("Invalid column selected for grouping.")
			}
		} else {
  		stop("Invalid sample grouping. Must be either a factor of length nrow(sc.counts)
  		     with two levels indicating training and validation set or name of a column in pheno that contains
			 entries 'training' or 'test' for each cell.
			")
		}
  	}
	}
	if (!all(is.logical(c(
	  simulation.bulks, simulation.genes, simulation.samples, simulation.subtypes
	 )))) {
		stop("Invalid value for at least one benchmark indicator.
		     Have to be logical.")
	}
	if (!is.numeric(repeats)) {
		stop("Invalid number of repeats. Must be numeric.")
	}
	if (!is.numeric(n.bulks)) {
		stop("Invalid input. n.bulks must be numeric.")
	}
	if (!is.numeric(n.cluster.sizes)) {
		stop("Invalid input. n.cluster.sizes must be integer")
	}
	if (repeats < 1) {
		warning("Repeats must be greater than 0. Setting to 1.")
		repeats <- 1
	}
	if (any(n.cluster.sizes < 1) | any(as.integer(n.cluster.sizes) != n.cluster.sizes)) {
		stop("n.cluster.sizes must be integer > 0")
	}
	if (n.bulks < 1) {
		warning("Number of bulks must be greater than 0. Setting to default.")
		n.bulks <- 1000
	}
	if (!is.null(exclude.from.bulks)) {
		if (!all(exclude.from.bulks %in% unique(sc.pheno[[cell.type.column]]))) {
			stop("Unknown cell type(s) in exclude.from.bulks.
			     Please select only cell types present in pheno data.")
		}
	}
	if (!is.null(exclude.from.signature)) {
		if (!all(exclude.from.signature %in% unique(sc.pheno[[cell.type.column]]))) {
			stop("Unknown cell types(s) in exclude.from.signature.
			     Please select only cell types present in pheno data.")
		}
	}

	# check and process algorithms input
	algorithms <- list(
			   list(algorithm = run_dtd, name = "DTD", model = NULL), 
			   list(
			     algorithm = run_least_squares, name = "Least_Squares", model = NULL
			   )
	)
	algorithm.names <- sapply(algorithms, function(x) {x$name})

	# create input algorithm list
	if (!is.null(input.algorithms)) {
		# input.algorithms must be a list
		if (!is.list(input.algorithms)) {
			stop("Invalid algorithm input. input.algorithms must be a list.")
		}
		predef.algos <- c()
		new.algos <- list()
		# allow either lists containing 2-3 entries (for function, name, model)
		# or characters specifying one of the predefined algorithms
		for (a in input.algorithms) {
			if (is.list(a)) {
				# for now check only whether algorithm exists, not its output
				if (is.function(a$algorithm) && is.character(a$name)) {
					new.algos <- c(new.algos, list(list(
            algorithm = a$algorithm, name = a$name, model = a$model
          )))
				}else{
					if (!is.character(a$name)) {
						stop("Invalid algorithm.")
					}else{
						stop(paste0("Invalid algorithm: ", a$name))
					}
				}
			}else{
				if (is.character(a) && a %in% algorithm.names) {
					predef.algos <- c(predef.algos, which(algorithm.names == a))
				}else{
					stop("Invalid algorithm.
					     User supplied algorithms can be given as a list with two entries:
					     name (name of the algorithm) and algorithm
					     (wrapper to call the algorithm)")
				}
			}
		}
		# join supplied (new) algorithms and specified existing algorithms
		if (length(predef.algos) > 0 || length(new.algos) > 0) {
			if (length(predef.algos) > 0) {
				algorithms <- algorithms[predef.algos]
				if (length(new.algos) > 0) {
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
	names(algorithms) <- algorithm.names
	
	check_algorithms(algorithms)
	
	# Data preparation
	cat("data preparation...\t\t", as.character(Sys.time()), "\n", sep = "")
	if (!is.null(sc.counts)) {
  	if (length(class(sc.counts)) > 1) {
  	  cat("Converting count matrix to sparse matrix...\n")
  		sc.counts <- Matrix(sc.counts, sparse = TRUE)
  		class(sc.counts) <- "dgCMatrix"
  	} else {
  	  if (class(sc.counts) != "dgCMatrix") {
  	    cat("Converting count matrix to sparse matrix...\n")
  	    sc.counts <- Matrix(sc.counts, sparse = TRUE)
  	    class(sc.counts) <- "dgCMatrix"
  	  }
  	}
	}
	if (!is.null(bulk.counts)) {
	  if (length(class(bulk.counts)) == 1) {
	    if (class(bulk.counts) != "dgCMatrix") {
	      bulk.counts <- Matrix(bulk.counts, sparse = TRUE)
	      class(bulk.counts) <- "dgCMatrix"
	    }
	  } else if (length(class(bulk.counts)) > 1) {
	    bulk.counts <- Matrix(bulk.counts, sparse = TRUE)
	    class(bulk.counts) <- "dgCMatrix"
	  }
		
	}

	# load / process / store data
	# if it exists load previously processed data from temp
	if (file.exists(paste(output.folder, "/input_data/training_set.h5", sep = ""))
      && file.exists(paste(output.folder, "/input_data/validation_set.h5",
         sep = ""))) {
		if (verbose) cat("Loading data found in temp directory\n")
		training_set <- read_data(
      paste(output.folder, "/input_data/training_set.h5", sep = "")
    )
		validation_set <- read_data(
      paste(output.folder, "/input_data/validation_set.h5", sep="")
    )
		training.exprs <- training_set$sc.counts
		training.pheno <- training_set$sc.pheno
		test.exprs <- validation_set$sc.counts
		test.pheno <- validation_set$sc.pheno
		sim.bulks <- list(
		  bulks = validation_set$bulk.counts,
		  props = validation_set$bulk.props
		)
		rm(list = c("training_set", "validation_set"))
	}else{
		dir.create(paste(output.folder, "/input_data", sep = ""), recursive = TRUE)
	}

	# save input data of benchmark() to temp directory
	function.call <- match.call()
	if (!file.exists(paste(output.folder, "input_data/raw.h5", sep = "/"))) {
	  suppressMessages(suppressWarnings(
      write_data(
        sc.counts, sc.pheno, bulk.counts, bulk.props,
        filename = paste(output.folder,"input_data/raw.h5", sep = "/")
      )
    ))
	}
	if (!file.exists(paste(output.folder, "input_data/params.h5", sep = "/"))) {
	  suppressMessages(suppressWarnings(
      write_misc_input(
        algorithm.names = algorithm.names, genesets = genesets,
        function.call = function.call, grouping = grouping,
        file = paste(output.folder, "input_data/params.h5", sep = "/")
      )
    ))
	}else{
	  # join list of algorithms and genesets from previous function calls
	  input_params <- read_misc_input(
      paste(output.folder, "input_data/params.h5", sep = "/")
    )
	  suppressMessages(suppressWarnings(
      write_misc_input(
        algorithm.names = unique(c(algorithm.names, input_params$algorithms)),
        genesets = unique(c(genesets, input_params$genesets)),
        function.call = function.call, grouping = grouping,
        file = paste(output.folder, "input_data/params.h5", sep = "/")
      )
    ))
	  rm("input_params")
	}

	# if any of the required data is missing preprocess input data
	if (!exists("training.exprs") || !exists("training.pheno") ||
      !exists("test.exprs") || !exists("test.pheno") || !exists("sim.bulks")) {
		if (verbose) cat("preprocessing data\n")
		if (cpm) {
			if(verbose) cat("scaling expression profiles to cpm\n")
      if (!is.null(sc.counts)) {
			     sc.counts <- scale_to_count(sc.counts)
      }
			if (!is.null(bulk.counts)) {
				bulk.counts <- scale_to_count(bulk.counts)
      }
      gc()
		}

		if (!is.null(grouping)) {
  		if(verbose) cat("splitting into test and training data\n")
  		# split data into test and validation set
  		if (length(unique(grouping)) == 2 && all(unique(grouping) %in% c(1, 2))) {
  			split.data <- split_dataset(
          sc.counts, sc.pheno, method = "predefined", grouping = grouping
        )
  			training.exprs <- split.data$training$exprs
  			training.pheno <- split.data$training$pheno
  			test.exprs <- split.data$test$exprs
  			test.pheno <- split.data$test$pheno

  			# remove unnecessary variable
  			rm("split.data")
  		}else if (length(unique(grouping)) == 1) {
  			message("Grouping vector contains only one group.
                 Disabling all simulations.")
  			training.exprs <- sc.counts
  			training.pheno <- sc.pheno
  			test.exprs <- NULL
  			test.pheno <- NULL
  			simulation.bulks <- FALSE
  			simulation.genes <- FALSE
  			simulation.samples <- FALSE
  			simulation.subtypes <- FALSE
  		}else{
  			stop(
  			  "Invalid grouping vector. Must be a vector of 
  			  1 (training samples) and 2 (test samples)."
  			  )
  		}
		}else{
			training.exprs <- NULL
			test.exprs <- NULL
			training.pheno <- NULL
			test.pheno <- NULL
			simulation.bulks <- FALSE
			simulation.genes <- FALSE
			simulation.samples <- FALSE
			simulation.subtypes <- FALSE
     		}

	  # this is not needed any more
	  rm(list = c("sc.counts", "sc.pheno", "grouping"))
	  gc()

	  include.in.bulks <- NULL
		# exclude.from.bulks is a parameter of benchmark(),
    # create_bulks expects the complement
		if (!is.null(exclude.from.bulks) &&
        length(intersect(
          unique(test.pheno[[cell.type.column]]),
          exclude.from.bulks
        )) > 0) {
			include.in.bulks <- unique(test.pheno[[cell.type.column]])[
        -which(unique(test.pheno[[cell.type.column]] %in% exclude.from.bulks))
      ]

		}
		# create simulated bulks if test data is available and they are needed;
	  # only needed if bulk/geneset/training set simulations are performed
	  # subtype simulation generates its own bulks
		if (!is.null(test.exprs) &&
        any(c(simulation.bulks, simulation.genes, simulation.samples))
      ) {
			cat("creating artificial bulks for simulation\n")
			sim.bulks <- create_bulks(
			  test.exprs,
			  test.pheno,
			  n.bulks,
			  include.in.bulks,
			  sum.to.count = cpm,
			  cell.type.column = cell.type.column,
			  n.profiles.per.bulk = n.profiles.per.bulk,
			  frac = 15
			)
		}else{
			sim.bulks <- list(bulks = NULL, props = NULL)
		}
		# save processed data for later use and repeated benchmarks
		suppressMessages(suppressWarnings(
      write_data(
        training.exprs, training.pheno,
        filename = paste(output.folder, "/input_data/training_set.h5", sep = "")
      )
    ))
		suppressMessages(suppressWarnings(
      write_data(
        test.exprs, test.pheno, sim.bulks$bulks, sim.bulks$props,
        filename = paste(
          output.folder, "/input_data/validation_set.h5", sep = ""
        )
      )
    ))
	}

	if(!is.null(training.pheno) && !is.null(training.exprs)){
  	# assume that samples in expression and pheno data are in the correct order
  	# if names do not match, assign sample names from expression to pheno data
  	if (!identical(rownames(training.pheno), colnames(training.exprs))) {
  		rownames(training.pheno) <- colnames(training.exprs)
  		colnames(training.exprs) <- rownames(training.pheno)
  	}
	}
	if (!is.null(test.exprs)) {
	  if (!identical(rownames(test.pheno), colnames(test.exprs))) {
	    rownames(test.pheno) <- colnames(test.exprs)
	    colnames(test.exprs) <- rownames(test.pheno)
	  }
	}

	if (!is.null(training.pheno)) {
  	# make sure subtypes of types in exclude.from.signature are also excluded
  	if (any(training.pheno[[cell.type.column]] %in% exclude.from.signature) &&
        "subtype" %in% colnames(training.pheno)) {
          subtype_names <- paste(
            training.pheno[[cell.type.column]],
            training.pheno$subtype,
            sep = "."
          )
          subtypes_exclude <- which(
            training.pheno[[cell.type.column]] %in% exclude.from.signature
          )
          exclude.from.signature <- c(
            exclude.from.signature,
            unique(subtype_names[subtypes_exclude])
          )
  	}

  	if (any(test.pheno[[cell.type.column]] %in% exclude.from.signature) &&
        "subtype" %in% colnames(test.pheno)) {
           subtype_names <- paste(
             test.pheno[[cell.type.column]] ,
             test.pheno$subtype,
             sep = "."
           )
           subtypes_exclude <- which(
             test.pheno[[cell.type.column]] %in% exclude.from.signature
           )
  	       exclude.from.signature <- c(
             exclude.from.signature,
             unique(subtype_names[subypes_exclude])
           )
  	}
  	exclude.from.signature <- unique(exclude.from.signature)
	}

	# deconvolution
	cat("deconvolution...\t\t", as.character(Sys.time()), "\n", sep = "")

	# select algorithms that have not been evaluated in a previous run
  result_dir <- paste(output.folder,"/results/real/",sep="")
   
  present.algorithms <- present_algos(
    target_dir = result_dir,
    name_pattern = "deconv.*.h5"
  )
  
  if (is.null(present.algorithms)) {
    to.run <- seq_len(length(algorithms))
  }else{
		if (verbose) cat("Found results for: ", present.algorithms, "\n")
		to.run <- which(! algorithm.names %in% present.algorithms)
	}

	res.no <- length(list.files(path = result_dir, pattern = "deconv.*.h5")) + 1
	if (!is.null(bulk.counts) && !is.null(bulk.props)) {
  	# deconvolute real bulks
  	cat("deconvolve real bulks...\t", as.character(Sys.time()), "\n", sep = "")
  	if (length(to.run) > 0) {
  		if (verbose) {
        cat(
          "algorithms to run:",
          sapply(algorithms[to.run], function(x) {x$name}),
          "\n", sep = " "
        )
      }
  		real.benchmark <- deconvolute(
  			training.expr = training.exprs,
  			training.pheno = training.pheno,
  			test.expr = NULL, test.pheno = NULL,
  			algorithms = algorithms[to.run],
  			verbose = verbose,
        exclude.from.bulks = NULL,
        exclude.from.signature = exclude.from.signature,
        max.genes = NULL,
        n.bulks = 0,
        bulks = list(bulks = bulk.counts, props = bulk.props),
        n.repeats = repeats,
        cell.type.column = cell.type.column,
        patient.column = patient.column
  		)
      # write deconvolution results
  		suppressMessages(suppressWarnings(
        write_result_list(
          real.benchmark,
          filename = paste(
            output.folder,
            "/results/real/deconv_output_",
            res.no, ".h5", sep = ""
          )
        )
      ))

  		# perform bootstrapping
  		cat("bootstrapping...\t\t", as.character(Sys.time()), "\n", sep = "")
  		estimates <- list()
  		# use the results of the algorithms in the first repetition
  		for (a in algorithms[to.run]) {
        tmp.result <- real.benchmark$results.list[["1"]][[a$name]]
  			estimates[[a$name]] <- tmp.result$est.props
  		}
  		props <- list(real = bulk.props, est = estimates)
  		bootstrap.real <- bootstrap_bulks(props)
  		h5_write_mat(
        bootstrap.real,
        filename = paste(
          output.folder,
          "/results/real/bootstrap_bulks",
          res.no, ".h5", sep = ""
        )
      )

      # remove real bulk results from RAM
    	rm(list = c("estimates", "props", "bootstrap.real"))
    	gc()
  	}
	}else{
	  real.benchmark <- NULL
	}

	# iterate through supplied simulation vector and perform those that are TRUE
	available.sims <- c(
    simulation.genes, simulation.samples, simulation.bulks, simulation.subtypes
  )
	names(available.sims) <- c("genes", "samples", "bulks", "subtypes")
	if (any(available.sims)) {
    cat("starting simulations...\t\t", as.character(Sys.time()), "\n", sep = "")
  }
	for (s in names(available.sims)) {
		if (!available.sims[s]) next
		# read previous results and exclude present algorithms

    sim_dir <- paste(output.folder, "/results/simulation/", s, "/", sep = "")
    if (s != "scores") {
      present.algorithms <- present_algos(sim_dir, name_pattern = "*.h5")
    }else{
      filename <- paste0(sim_dir, "deconv_output.h5")
      if (file.exists(filename)) {
        prev_res <- rhdf5::h5dump(filename)
        prev_runs <- length(prev_res)
        res.no <- prev_runs + 1
        if (prev_runs > 0) {
          present.algorithms <- unique(unlist(
            sapply(
              prev_res,
              FUN = function(x) {
                unique(x$algorithm)
              }
            )
          ))
        }
      }else{
        res.no <- 0
        present.algorithms <- NULL
      }
    }
    
    if (is.null(present.algorithms)) {
      to.run <- seq_len(length(algorithms))
    }else{
			if(verbose) print(present.algorithms)
			to.run <- which(! algorithm.names %in% present.algorithms)
		}
    rm("present.algorithms")
    if (s != "scores") {
		  res.no <- length(list.files(path = sim_dir, pattern = "*.h5")) + 1
    }

		# execute benchmark corresponding to s and save results
		if (length(to.run) > 0) {
			if (s == "bulks") {
				cat("bulk simulation...\t\t", as.character(Sys.time()), "\n", sep = "")
				benchmark.results <- deconvolute(
          training.expr = training.exprs,
          training.pheno = training.pheno,
          test.expr = NULL, test.pheno = NULL,
          algorithms = algorithms[to.run],
          verbose = verbose,
          exclude.from.bulks = NULL,
          exclude.from.signature = exclude.from.signature,
          max.genes = NULL, n.bulks = 0,
          bulks = sim.bulks,
          n.repeats = repeats,
          cell.type.column = cell.type.column,
          patient.column = patient.column
        )
			}
			if (s == "genes") {
				cat(
          "geneset simulation...\t\t", as.character(Sys.time()), "\n", sep = ""
        )
			  benchmark.results <- geneset_benchmark(
          training.exprs = training.exprs,
          training.pheno = training.pheno,
          test.exprs = NULL, test.pheno = NULL,
          genesets = genesets,
          algorithms = algorithms[to.run],
          bulk.data = sim.bulks,
          n.repeats = repeats,
          exclude.from.signature = exclude.from.signature,
          verbose = verbose,
          cell.type.column = cell.type.column,
          patient.column = patient.column
        )
			}
			if (s == "samples") {
				cat(
          "sample simulation...\t\t", as.character(Sys.time()), "\n", sep = ""
        )
			  benchmark.results <- sample_size_benchmark(
          training.exprs = training.exprs,
          training.pheno = training.pheno,
          test.exprs = NULL,
          test.pheno = NULL,
          algorithms = algorithms[to.run],
          bulk.data = sim.bulks,
          n.repeats = repeats,
          exclude.from.signature = exclude.from.signature,
          step.size = 0.25,
          verbose = verbose,
          cell.type.column = cell.type.column,
          patient.column = patient.column
        )
			}
			if (s == "subtypes") {
				cat(
          "subtype simulation...\t\t", as.character(Sys.time()), "\n", sep = ""
        )
				if (!is.null(test.exprs) && !is.null(test.pheno)) {
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
					n.clusters = c(1, 2, 4, 8),
					patient.column = patient.column,
					n.bulks = n.bulks,
          repeats = repeats
				)

			  rm(list = c("all.exprs", "all.pheno"))
			  gc()
			}

			if (!dir.exists(
            paste(output.folder, "/results/simulation/", s, sep = "")
      )) {
		 		dir.create(
          paste(output.folder, "/results/simulation/", s, sep = ""),
          recursive = TRUE
        )
			}
        suppressMessages(suppressWarnings(
          write_result_list(
            benchmark.results,
            filename = paste(
              output.folder,
              "/results/simulation/",
              s, "/deconv_output_", res.no, ".h5", sep = ""
            )
          )
        ))
      
      rm("benchmark.results")
    }
	}
	if (verbose) cat(
	  "Creating plots...\t\t", as.character(Sys.time()), "\n", sep = ""
	  )
	celltype_order_list <- plot_all(
	  temp_dir = output.folder,
	  genesets = genesets,
	  features = rownames(training.exprs)
	)
	
	if (report) {
  	cat("preparing results...\t\t", as.character(Sys.time()), "\n", sep = "")
  	# deconvolution step is over
  	# create results markdown
  	report.path <- suppressWarnings(
      render_results(
        output.folder, benchmark.name, cell.type.column, 
		celltype_order_list$celltype_order, celltype_order_list$celltype_order_simulated
      )
    )
  	cat("Done\t\t\t\t", as.character(Sys.time()), "\n\n", sep = "")
  	cat("Report generated: ", report.path, "\n", sep = "")
  	cat("Created plots can be found in: ", output.folder, "/report_plots/\n", sep = "")
	}else{
	  report.path <- NULL
	}
	if(dir.exists("CIBERSORT")){
	  unlink("CIBERSORT", recursive = TRUE)
	}

	if (verbose) {
		t <- tictoc::toc(quiet = TRUE)
		time <- t$toc - t$tic
		if(time > 3600) {
			time <- time / 3600
			cat(
        t$msg, ": ", as.integer(time), "h ",
        as.integer((time-as.integer(time))*60), "m elapsed\n", sep = ""
      )
		}else{
			if(time > 60) {
				time <- time / 60
				cat(t$msg, ": ", as.integer(time), "m ",
            as.integer((time-as.integer(time))*60), "s elapsed\n", sep = ""
        )
			} else {
				cat(t$msg, ": ", time, "s elapsed\n", sep = "")
			}
		}
	}
	if (exists("real.benchmark")) {
		return(list(report_path = report.path, bulk_results = real.benchmark))
	} else {
		return(list(report_path = report.path, bulk_results = NULL))
	}
}
