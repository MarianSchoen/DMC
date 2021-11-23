#' create condition number plots for deconvolution results
#'
#' @param results.df data frame containing results as returned by prepare_data
#' @param algorithm.order character vector of algorithm names
#' specifying the plotting order
#' @return list containing 3 plots:\cr
#' 1) cond_num_plot - barplot of average condition numbers\cr
#' 2) cond_vs_score - scatter plot of
#' average condition number vs. score for each algorithm\cr
#' 3) variation_plot - scatter plot of standard deviation of
#' score vs. standard deviation of condition number

plot_cond_num <- function(
  results.df,
  algorithm.order = NULL
) {
  if (!is.data.frame(results.df)) {
      stop("results.df must be a data frame")
  }
  required_cols <- c("score", "algorithm", "cell_type", "condition_number")
  if (!all(required_cols %in% colnames(results.df))) {
      stop("required columns missing from results.df")
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
  # use only 'overall' rows
  results.df <- results.df[which(results.df$cell_type == "overall"), ]
  overall.df <- data.frame()

  # create means and sds of scores and condition numbers from repetitions
  for (a in unique(results.df$algorithm)) {
      overall.df <- rbind(
        overall.df,
        data.frame(
          algorithm = a,
          condition_number = mean(
            results.df[which(results.df$algorithm == a), "condition_number"]
          ),
          condition_variation = sd(
            results.df[which(results.df$algorithm == a), "condition_number"]
          ),
          score = mean(results.df[which(results.df$algorithm == a), "score"]),
          score_variation = sd(
            results.df[which(results.df$algorithm == a), "score"]
          )
        )
      )
  }

  # order the factor levels for plot order
  if (!is.null(algorithm.order)) {
    overall.df$algorithm <- factor(
      overall.df$algorithm,
      levels = algorithm.order
    )
  }
  # create plot labels (for barplot)
  cond_labs <- levels(overall.df$algorithm)
  for (i in seq_len(length(cond_labs))) {
      idx <- which(overall.df$algorithm == cond_labs[i])
      if (is.na(overall.df[idx, "condition_number"])) {
          cond_labs[i] <- paste(cond_labs[i], "\n data missing", sep = "")
      }
  }

  # plot condition numbers as barplot
  cond_num_plot <- ggplot(overall.df) +
    geom_bar(
      aes(x = algorithm, y = condition_number, fill = algorithm),
      stat = "identity",
      position = "dodge"
    ) +
    geom_errorbar(
      aes(
        x = algorithm,
        ymin = condition_number - condition_variation,
        ymax = condition_number + condition_variation
      ),
      width = 0.2
    ) +
    ggtitle("average signature matrix condition number") +
    ylab("condition number") +
    xlab("algorithm") +
    scale_x_discrete(
      limits = levels(overall.df$algorithm), labels = cond_labs
    ) +
    theme(
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        title = element_text(size = 24),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20)
    )

    # plot score vs condition number
  cond_vs_score <- ggplot(overall.df) +
    geom_point(
      aes(x = condition_number, y = score, col = algorithm),
      size = 3, na.rm = TRUE
    ) +
    ggtitle("score vs condition number") +
    xlab("condition number") +
    ylab("score") + ylim(0, 1) +
    theme(
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        title = element_text(size = 24),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20)
    )

  # plot sd of score vs sd of condition number
  variation_plot <- ggplot(overall.df) +
    geom_point(
      aes(x = condition_variation, y = score_variation, col = algorithm),
      size = 3, na.rm = TRUE
    ) +
    ggtitle("SDs of score vs SDs of condition number") +
    xlab("condition number SD") +
    ylab("score SD") +
    theme(
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        title = element_text(size = 24),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20)
    )

  return(list(
    cond_num_plot = cond_num_plot,
    cond_vs_score = cond_vs_score,
    variation_plot = variation_plot
  ))
}
