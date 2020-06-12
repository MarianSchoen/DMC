#' render markdown report for data produced by \link{benchmark}
#' 
#' @param temp.dir directory containing the benchmark data. \cr
#' Usually \code{paste(temp.dir,"/",benchmark.name,sep="")}, 
#' where \code{temp.dir} and \code{benchmark.name} 
#' are parameters of \link{benchmark}.
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @return NULL
#' @export

render_results <- function(temp.dir, metric = "cor"){
	# parameter checks
	if(!dir.exists(temp.dir)){
		stop("Invalid temp directory")
	}
	if(is.character(metric)){
    if(metric != "cor"){
			stop("metric must be either \"cor\" or a function")
		}else{
			metric <- cor
		}
  }else{
    if(!is.function(metric)){
			stop("Function corresponding to 'metric' could not be found.")
		}
  }
	
	# render the template to pdf with the data stored in temp.dir
	render(
	  input = system.file(
	    "rmd", 
	    "report.Rmd"
	    , package = "DeconvolutionAlgorithmBenchmarking")
	  , params = list(tempdir = temp.dir, metric=metric)
	  , output_file = paste(temp.dir, "/report_", gsub(" ", "_", Sys.time()),".pdf", sep = "")
	  )
}