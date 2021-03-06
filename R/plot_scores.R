#' create score plot
#'
#' @param results.df data frame as returned by prepare_data
#' @param title character, plot title
#' @param real.props non-negative numeric matrix, with cell types as rows,
#' and bulk RNA-Seq profiles.
#' @param celltype.order character vector of cell types specifying
#' the plotting order
#' @param algorithm.order character vector of algorithm names
#' specifying the plotting order (left to right)
#' @result list with 3 entries:
#' 1) plot - table plot of deconvolution results
#' 2) celltype.order - character vector containing cell types in plotting order
#' 3) algorithm.order - character vector containing algorithms in plotting order

evaluation_plot <- function(
  results.df,
  title = NULL,
  real.props = NULL,
  celltype.order = NULL,
  algorithm.order = NULL
) {
  # parameter check
  if (!is.data.frame(results.df)) {
    stop("results.df must be a data frame")
  }
  required_names <- c("algorithm", "score", "cell_type")
  if (!all(required_names %in% colnames(results.df))) {
    stop("required columns missing from results.df")
  }
  if (is.null(title)) {
      title <- "deconvolution quality"
  }
  
  if (!is.null(real.props) && !is.matrix(real.props)) {
    stop("real.props is not a matrix")
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
            stop("algorithm.order does not fit the
                algorithm column of results.df")
      }
  }
  
  # cannot plot negative values -> cutoff at 0
  results.df[which(!is.na(results.df$score) & results.df$score < 0), "score"] <- 0

  # reduce to one entry per cell type and algorithm
  # by taking the mean over repetitions
  quality.scores <- c()
  for (a in unique(results.df$algorithm)) {
    for (ct in unique(results.df$cell_type)) {
      subset <- which(results.df$algorithm == a & results.df$cell_type == ct)
      if (length(subset > 0)) {
        quality.scores <- rbind(
          quality.scores,
          c(
            mean(results.df[subset, ]$score),
            mean(results.df[subset, ]$error),
            sd(results.df[subset, ]$score),
            ct, a, "correlation"
          )
        )
      }
    }
  }

  if (is.null(nrow(quality.scores))) {
      warning("No scores to plot")
      return(NULL)
  }
  # convert to data frame for plotting
  quality.scores <- as.data.frame(quality.scores)
  colnames(quality.scores) <- c(
    "value", "bootstrap_error","sd", "cell_type", "algorithm", "measure"
  )
  # create labels for cell types in table plots
  if (is.null(celltype.order)) {
    labels <- unique(as.character(quality.scores$cell_type))
    # make sure 'overall' is the first type
    labels <- c("overall", labels[-which(labels == "overall")])
  }else{
    labels <- celltype.order
  }
  quality.scores$cell_type <- factor(quality.scores$cell_type, levels = labels)

  # if real.props is available,
  # add average fraction of this cell type in the bulks to cell type label
  if (!is.null(real.props)) {
    for (i in seq_len(length(labels))) {
      if (labels[i] != "overall") {
        prop <- round(mean(real.props[labels[i], ]), 2)
        labels[i] <- paste(labels[i], "\n", prop * 100, "%", sep = "")
      }
    }
  }

  
  # cleaner plot by setting negative correlations to 0
  quality.scores$value <- as.numeric(as.character(quality.scores$value))
  quality.scores$sd <- as.numeric(as.character(quality.scores$sd))
  non_NA_indices <- which(!is.na(quality.scores$value))
  if (any(quality.scores$value[non_NA_indices] < 0)) {
    quality.scores$value[non_NA_indices[non_NA_indices < 0]] <- 0
  }

  # order algorithms by performance if ordering is not given
  algos <- quality.scores$algorithm[
    which(quality.scores$cell_type == "overall")
  ]
  ranking <- quality.scores$value[which(quality.scores$cell_type == "overall")]
  sds <- quality.scores$sd[which(quality.scores$cell_type == "overall")]
  names(ranking) <- algos
  names(sds) <- algos

  if (is.null(algorithm.order)) {
    quality.scores$algorithm <- factor(
      quality.scores$algorithm,
      levels = names(ranking)[order(ranking)]
    )
    sds <- sds[order(ranking)]
    ranking <- sort(ranking)
  }else{
    # order by algorithm.order
    quality.scores$algorithm <- factor(
      quality.scores$algorithm,
      levels = algorithm.order
    )
    sds <- sds[algorithm.order]
    ranking <- ranking[algorithm.order]
  }

  # create table plot
  score.plot <- ggplot(
    quality.scores,
    aes(x = as.numeric(cell_type), y = as.numeric(algorithm))
  ) +
  geom_point(aes(size = value, col = value, fill = value), shape = 22)  +
  geom_point(shape = 22, fill = NA, color = "black", size = 25) +
  xlab("cell type") + ylab("algorithm") +
  ggtitle(title, subtitle = "correlation") +
  scale_size_continuous(limits = c(0, 1), range = c(0, 25), guide = "none") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(
    breaks = seq_len(length(levels(quality.scores$algorithm))),
    limits = c(0.5, length(levels(quality.scores$algorithm)) + 0.5),
    labels = paste(
      levels(quality.scores$algorithm),
      "\nscore=",
      round(ranking, 2),
      "\n+/- ",
      round(sds, 2),
      sep = ""
    ),
    minor_breaks = seq(
      0., length(levels(quality.scores$algorithm)) + 1, 0.01
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(length(levels(quality.scores$cell_type))),
    limits = c(0.5, length(levels(quality.scores$cell_type)) + 0.5),
    labels = labels,
    minor_breaks = seq(
      0., length(levels(quality.scores$cell_type)) + 1., 0.01
    ),
    sec.axis = dup_axis()
  ) +
  labs(fill = "Score") +
  scale_color_gradient2(
    low = "red", mid = "orange", high = "green",
    midpoint = 0.5, limits = c(0, 1), guide = "none"
  ) +
  scale_fill_gradient2(
    low = "red", mid = "orange", high = "green",
    midpoint = 0.5, limits = c(0,1)
  ) +
  geom_vline(xintercept = 1.5) +
  coord_cartesian(
    xlim = c(0.75, length(levels(quality.scores$cell_type)) + 0.25)
  )
  # return the plot and the ordering of the axes
  return(list(
    plot = score.plot,
    celltype.order = levels(quality.scores$cell_type),
    algorithm.order = levels(quality.scores$algorithm)
  ))
}
