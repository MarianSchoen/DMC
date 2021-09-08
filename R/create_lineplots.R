#' create lineplots of deconvolution results and runtime for different gene sets
#'
#' @param results.df data frame as returned by prepare_data
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @param metric.name string, name of the evaluation metric used;
#' not needed if metric is a string ("cor").
#' If metric is a function and metric.name
#' is NULL, the default will be "custom metric"
#' @param genesets list of gene sets (character vectors)
#' @param available.features character vector containing names
#' of available features
#' @param celltype.order character vector of cell types
#' specifying the plotting order
#' @param algorithm.order character vector of algorithm names
#' specifying the plotting order
#' @return list containing two elements:\cr
#' 1) runtime.plot\cr
#' 2) list containing lineplot of scores for each celltypes\cr

create_lineplots <- function(
  results.df,
  metric = "cor",
  metric.name = NULL,
  genesets = NULL,
  available.features = NULL,
  celltype.order = NULL,
  algorithm.order = NULL
) {
  # parameter checks
  if (!is.data.frame(results.df)) {
    stop("results.df must be a data frame")
  }
  required_cols <- c(
    "algorithm", "score", "metric", "geneset", "cell_type", "time"
  )
  if (!all(required_cols %in% colnames(results.df))) {
    stop("required columns missing from results.df")
  }
  if (!is.null(genesets)) {
    if (!is.list(genesets)) {
        stop("genesets must be a named list of character vectors")
    }
  }
  if (!is.null(available.features)) {
    if (!is.character(available.features)) {
      stop("available features must be a charcter vector")
    }
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
    if (!is.function(metric)) {
      stop("Function corresponding to 'metric' could not be found.")
    }else{
      if (is.null(metric.name) || !is.character(metric.name)) {
        metric.name <- "custom metric"
      }
    }
  }

  overall.df <- results.df[which(results.df$cell_type == "overall"), ]
  # order algorithms by performance or given order
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
  if (!is.null(celltype.order)) {
    results.df$cell_type <- factor(
      results.df$cell_type,
      levels = celltype.order
    )
  }


  # create display labels for gene sets and
  # sort according to number of genes if possible
  if (!is.null(genesets)) {
    if (all(unique(results.df$geneset) %in% names(genesets))) {
      geneset.labs <- paste(
        names(genesets),
        "\n(",
        as.numeric(
          sapply(genesets, function(x) length(which(x %in% available.features)))
        ),
        " genes)",
        sep = ""
      )
    }else{
      geneset.labs <- levels(results.df$geneset)
    }
    geneset.sizes <- sapply(
      genesets,
      function(x) length(which(x %in% available.features))
    )
    geneset.labs <- geneset.labs[order(geneset.sizes)]
    geneset.limits <- names(sort(geneset.sizes))
  }else{
    geneset.labs <- levels(results.df$geneset)
    geneset.limits <- levels(results.df$geneset)
  }

  # create plot per cell type (including overall)
  cell.type.plots <- list()
  for (t in levels(results.df$cell_type)) {
    sub.df <- results.df[which(results.df$cell_type == t), ]
    # create data frame containing score for each
    # algorithm and gene set for this cell type
    temp.scores <- tapply(
      sub.df$score,
      list(sub.df$algorithm, sub.df$geneset),
      mean
    )
    temp.sds <- tapply(
      sub.df$score,
      list(sub.df$algorithm, sub.df$geneset),
      sd
    )
    temp.times <- tapply(
      sub.df$time,
      list(sub.df$algorithm, sub.df$geneset),
      mean
    )
    temp.df <- data.frame()
    for (i in seq_len(ncol(temp.scores))) {
      for (j in seq_len(nrow(temp.scores))) {
        temp.df <- rbind(
          temp.df,
          data.frame(
            algorithm = rownames(temp.scores)[j],
            geneset = colnames(temp.scores)[i],
            score = temp.scores[j, i],
            time = temp.times[j,i],
            sd = temp.sds[j, i]
          )
        )
      }
    }
    temp.df$geneset <- factor(temp.df$geneset, levels = geneset.limits)

    cell.type.plots[[t]] <- ggplot(
      temp.df,
      aes(x=geneset, y = score, group = algorithm, col = algorithm)
    ) +
    geom_line(size = 2) +
    geom_point() +
    geom_errorbar(
      aes(x = geneset, ymin = score - sd, ymax = score + sd),
      width = 0.2
    ) +
    xlab("gene set (increasing size)") +
    ylab(metric.name) +
    ggtitle(paste(
      "deconvolution quality using different gene sets (", t, ")", sep = ""
    )) +
    scale_x_discrete(limits = geneset.limits, labels = geneset.labs) +
    guides(linetype = guide_legend(override.aes = list(size = 2))) +
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22),
      title = element_text(size = 24),
      axis.title.x = element_text(size = 22),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 20)
    )
    cell.type.plots[[t]] <- cell.type.plots[[t]] + ylim(0, 1)
    # create only one runtime plot
    if (t == "overall") {
      runtime.plot <- ggplot(
        temp.df,
        aes(x = geneset, y = log(time, 10), group = algorithm, col = algorithm)
      ) +
      geom_line(size = 2) +
      geom_point() +
      xlab("gene set") +
      ylab("log time (s)") +
      ggtitle("runtime of algorithms using different gene sets") +
      theme(
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        title = element_text(size = 24),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20)
      ) +
      scale_x_discrete(limits = geneset.limits, labels = geneset.labs) +
      guides(linetype = guide_legend(override.aes = list(size = 2)))
    }
  }
  return(list(
    runtime.plot = runtime.plot,
    cell.type.plots = cell.type.plots
  ))
}
