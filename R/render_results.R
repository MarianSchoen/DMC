render_results <- function(temp.dir, metric){
	require(rmarkdown)
	render(
	  input = system.file(
	    "rmd", 
	    "report.Rmd"
	    , package = "DeconvolutionAlgorithmBenchmarking")
	  , params = list(tempdir = temp.dir, metric=metric)
	  , output_file = paste(temp.dir, "/report_",Sys.time(),".pdf", sep = "")
	  )
}