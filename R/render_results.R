#' render markdown report for data produced by \link{benchmark}
#' 
#' @param temp.dir directory containing the benchmark data. \cr
#' Usually \code{paste(temp.dir,"/",benchmark.name,sep="")}, 
#' where \code{temp.dir} and \code{benchmark.name} 
#' are parameters of \link{benchmark}.
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @param metric.name string, name of the evaluation metric used; not needed if metric is a string ("cor"). If metric is a function and metric.name
#' is NULL, the default will be "custom metric"
#' @return NULL
#' @export

render_results <- function(temp.dir, metric = "cor", metric.name = NULL){
	# parameter checks
	if(!dir.exists(temp.dir)){
		stop("Invalid temp directory")
	}
	if(is.character(metric)){
    if(metric != "cor"){
			stop("metric must be either \"cor\" or a function")
		}else{
			if(is.null(metric.name) || !is.character(metric.name)){
				metric.name <- "custom metric"
			}
			metric <- cor
		}
  }else{
    if(!is.function(metric)){
			stop("Function corresponding to 'metric' could not be found.")
		}else{
			if(is.null(metric.name) || !is.character(metric.name)){
				metric.name <- "custom metric"
			}
		}
  }
	
	# render the template to pdf with the data stored in temp.dir
	render(
	  input = system.file(
	    "rmd", 
	    "report.Rmd"
	    , package = "DeconvolutionAlgorithmBenchmarking")
	  , params = list(tempdir = temp.dir, metric=metric, metric.name = metric.name)
	  , output_file = paste(temp.dir, "/report_", gsub(" ", "_", Sys.time()),".pdf", sep = "")
	  )
}