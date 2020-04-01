render_results <- function(temp.dir, metric = "cor"){
	if(!dir.exists(temp.dir)){
		stop("Invalid temp directory")
	}
	if(!metric %in% c("cor", "mad", "rmsd")){
		stop("Invalid metric. Must be one of 'cor', 'mad', 'rmsd'")
	}
	library(rmarkdown)
	render(
	  input = system.file(
	    "rmd", 
	    "report.Rmd"
	    , package = "DeconvolutionAlgorithmBenchmarking")
	  , params = list(tempdir = temp.dir, metric=metric)
	  , output_file = paste(temp.dir, "/report_",Sys.time(),".pdf", sep = "")
	  )
}