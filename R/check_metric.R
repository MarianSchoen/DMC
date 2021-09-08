#' check and process the supplied metric and metric name
#' 
#' @param metric either "cor" or a function evaluating 
#' the similarity of two vectors
#' @param metric.name character containing name of the metric
#' @return list with entries metric, metric.name

check_metric <- function(metric, metric.name) {
  if (is.character(metric)) {
    if (metric != "cor") {
      stop("metric must be either \"cor\" or a function")
    }else{
      metric <- cor
      if (is.null(metric.name) || !is.character(metric.name)) {
        metric.name <- "pearson correlation"
      }
    }
  }else{
    if (is.function(metric)) {
      warning("Using custom evaluation function / metric.
			        Unexpected results and plots may occur.")
      if (is.null(metric.name) || !is.character(metric.name)) {
        metric.name <- "custom metric"
      }
    }else{
      stop("Function corresponding to 'metric' could not be found.")
    }
  }
  test.result <- metric(c(1,2,3), c(1,2,3))
  if (!is.numeric(test.result) || length(test.result) != 1) {
    stop("Invalid function. 'metric' must accept two vectors 
         and return one number (numeric) as result.")
  }
  return(list(metric = metric, metric.name = metric.name))
}