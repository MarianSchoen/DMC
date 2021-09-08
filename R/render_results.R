#' render markdown report for data produced by \link{benchmark}
#'
#' @param temp.dir directory containing the benchmark data. \cr
#' Usually \code{paste(temp.dir,"/",benchmark.name,sep="")},
#' where \code{temp.dir} and \code{benchmark.name}
#' are parameters of \link{benchmark}.
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @param metric.name string, name of the evaluation metric used;
#' not needed if metric is a string ("cor"). If metric is a function and
#' metric.name is NULL, the default will be "custom metric"
#' @param benchmark.name string, name of the benchmark to be included in
#' report filename
#' @param celltype.col string, name of the column in sc.pheno that contains
#' cell type information; default "cell_type"
#' @return NULL
#' @export

render_results <- function(
	temp.dir,
	metric = "cor",
	metric.name = NULL,
	benchmark.name = "",
	celltype.col = "cell_type"
) {
	# parameter checks
	if (!dir.exists(temp.dir)) {
		stop("Invalid temp directory")
	}
  # check metric
  metric.list <- check_metric(metric, metric.name)
  metric <- metric.list$metric
  metric.name <- metric.list$metric.name

	# render the template to pdf with the data stored in temp.dir
	rmarkdown::render(
	  input = system.file(
	   	"rmd",
	    "report.Rmd",
	    package = "DAB"
		),
  	params = list(
			tempdir = temp.dir,
			metric = metric,
			metric.name = metric.name,
			celltype.col = celltype.col
		),
	  output_file = paste(
			temp.dir, "/", benchmark.name, "_report_",
			gsub(" ", "_", Sys.time()), ".html", sep = ""
		),
	  quiet = TRUE
	)
}
