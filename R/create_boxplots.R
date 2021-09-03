#' create boxplots of deconvolution results for each algorithm
#'
#' @param results.df data frame as returned by prepare_data
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @param metric.name string, name of the evaluation metric used;
#' not needed if metric is a string ("cor").
#' If metric is a function and metric.name
#' is NULL, the default will be "custom metric"
#' @param celltype.order character vector of cell types
#' specifying the plotting order
#' @param algorithm.order character vector of algorithm names
#' specifying the plotting order (left to right)
#' @return list containing list of ggplot objects (plots),
#' ordered by celltype.order (if supplied)

create_boxplots <- function(
  results.df,
  metric = "cor",
  metric.name = NULL,
  celltype.order = NULL,
  algorithm.order = NULL
) {
  # parameter checks
  if (!is.data.frame(results.df)) {
      stop("results.df must be a data frame")
  }
  if (!all(c("score", "algorithm", "cell_type") %in% colnames(results.df))) {
      stop("required columns missing from results.df")
  }
  if (!is.null(celltype.order)) {
    if (!is.character(celltype.order)) {
        stop("celltype.order must be a charcter vector")
    }
    if (!all(celltype.order %in% unique(results.df$cell_type)) ||
        length(celltype.order) != length(unique(results.df$cell_type))) {
        stop("celltype.order does not fit the cell_type column of results.df")
    }
  }
  if (!is.null(algorithm.order)) {
    if (!is.character(algorithm.order)) {
        stop("celltype.order must be a charcter vector")
    }
    if (!all(algorithm.order %in% unique(results.df$algorithm)) ||
        length(algorithm.order) != length(unique(results.df$algorithm))) {
        stop("algorithm.order does not fit the algorithm column of results.df")
    }
  }
  if (is.character(metric)) {
    if (metric != "cor") {
      stop("metric must be either \"cor\" or a function")
    }else{
      if (is.null(metric.name) || !is.character(metric.name)) {
        metric.name <- "custom metric"
      }
      metric <- cor
    }
  }else{
    if (!exists(as.character(substitute(metric)))) {
      stop("Function corresponding to 'metric' could not be found.")
    }else{
      if (is.null(metric.name) || !is.character(metric.name)) {
        metric.name <- "custom metric"
      }
    }
  }
  overall.df <- results.df[which(results.df$cell_type == "overall"), ]

  # order algorithms by performance or by given vector
  if (is.null(algorithm.order)) {
    performances <- tapply(overall.df$score, overall.df$algorithm, median)
    results.df$algorithm <- factor(
      results.df$algorithm,
      levels = levels(overall.df$algorithm)[order(performances)]
    )
  }else{
    results.df$algorithm <- factor(
      results.df$algorithm,
      levels = algorithm.order
    )
  }
  # order cell types by supplied vector
  if (!is.null(celltype.order)) {
    results.df$cell_type <- factor(
      results.df$cell_type,
      levels = celltype.order
    )
  }

  # create boxplots for each cell type (including overall)
  cell.type.plots <- list()
  for (t in levels(results.df$cell_type)) {
    sub.df <- results.df[which(results.df$cell_type == t), ]
    cell.type.plots[[t]] <- ggplot(sub.df, aes(x = algorithm, y = score)) +
      geom_boxplot(aes(col = algorithm)) +
      xlab("algorithm") +
      ylab(metric.name) +
      ggtitle("quality of deconvolution results", subtitle = t) +
      theme(
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        title = element_text(size = 24),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20)
      ) +
      scale_x_discrete(limits = levels(results.df$algorithm))
    cell.type.plots[[t]] <- cell.type.plots[[t]] + ylim(0, 1)
    if (t == "overall") {
        cell.type.plots[[t]] <- cell.type.plots[[t]] +
          ylab(paste("average overall score"))
    }
  }
  return(list(cell.type.plots = cell.type.plots))
}
