#' render markdown report for data produced by \link{benchmark}
#' 
#' @param temp.dir directory containing the benchmark data. \cr
#' Usually \code{paste(temp.dir,"/",benchmark.name,sep="")}, 
#' where \code{temp.dir} and \code{benchmark.name} 
#' are parameters of \link{benchmark}.
#' @param metric character specifying the metric. default: 'cor'
#' @return NULL
#' @export

render_results <- function(temp.dir, metric = "cor"){
	# parameter checks
	if(!dir.exists(temp.dir)){
		stop("Invalid temp directory")
	}
	if(!metric %in% c("cor", "mad", "rmsd")){
		stop("Invalid metric. Must be one of 'cor', 'mad', 'rmsd'")
	}
	
	# render the template to pdf with the data stored in temp.dir
	render(
	  input = system.file(
	    "rmd", 
	    "report.Rmd"
	    , package = "DeconvolutionAlgorithmBenchmarking")
	  , params = list(tempdir = temp.dir, metric=metric)
	  , output_file = paste(temp.dir, "/report_",Sys.time(),".pdf", sep = "")
	  )
}