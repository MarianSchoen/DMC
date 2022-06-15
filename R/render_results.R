#' render markdown report for data produced by \link{benchmark}
#'
#' @param temp.dir directory containing the benchmark data. \cr
#' Usually \code{paste(temp.dir,"/",benchmark.name,sep="")},
#' where \code{temp.dir} and \code{benchmark.name}
#' are parameters of \link{benchmark}.
#' @param benchmark.name string, name of the benchmark to be included in
#' report filename
#' @param celltype.col string, name of the column in sc.pheno that contains
#' cell type information; default "cell_type"
#' @return NULL
#' @export

render_results <- function(
	temp.dir,
	benchmark.name = "",
	celltype.col = "cell_type"
) {
	# parameter checks
	if (!dir.exists(temp.dir)) {
		stop("Invalid temp directory")
	}

	# render the template to pdf with the data stored in temp.dir
	rmarkdown::render(
	  input = system.file(
	   	"rmd",
	    "report.Rmd",
	    package = "DMC"
		),
  	params = list(
			tempdir = temp.dir,
			celltype.col = celltype.col
		),
	  output_file = paste(
			temp.dir, "/", benchmark.name, "_report_",
			gsub(" ", "_", Sys.time()), ".html", sep = ""
		),
	  quiet = TRUE
	)
}
