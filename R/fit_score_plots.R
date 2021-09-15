#' This function creates fit plots and plots of offsets and errors (scores)
#' from the output of `score_algorithms`
#' 
#' @param score.results list of data frames as read from h5 result file; 
#' each data frame is one output of `score_algorithms`
#' @return list containing score plots:\cr
#' 1) fits: list containing one fit plot for each repetition\cr
#' 2) offset_plot: boxplot of fitted offset for each algorithm and celltype\cr
#' 3) error_plot: boxplot of fitted error for each algorithm and celltype\cr
#' 4) combined_plot: barplot of (offset + error) for each algorithm and celltype

create_fit_plots <- function (score.results) {
  if (!is.list(score.results) || length(score.results) == 0) {
    stop("In function fit_plots: Invalid input. Must be list of data frames.")
  }
  results.df <- score.results[[1]]
  if (length(score.results) > 1) {
    for (i in 2:length(score.results)) {
      results.df <- rbind(
        results.df,
        score.results[[i]]
      )
    }
  }

  # create fit plots
  fit_plots <- list()
  for (rep in seq_len(max(results.df$repetition))) {
    temp.df <- results.df[which(results.df$repetition == rep), ]
    fit_plots[[rep]] <- ggplot(temp.df, aes(x = variance)) +
      geom_point(aes(y = correlation)) +
      geom_line(aes(y = prediction), col = "red") +
      facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
      ggtitle(
        "Fitted correlation curves", 
        subtitle = paste0("repetition ", rep)
      ) +
      xlab("Variance in cellular composition") +
      ylab("Correlation")
  }

  # create plot of offsets (boxplot)
  offset_plots <- ggplot(results.df, aes(x = celltype, y = offset)) +
    geom_boxplot(aes(fill = celltype)) + 
    facet_grid(rows = vars(algorithm)) +
    ggtitle("fitted offset parameters") +
    xlab("cell type") +
    ylab("offset")

  # create plots of errors (boxplot)
  error_plots <- ggplot(results.df, aes(x = celltype, y = error_sd)) +
    geom_boxplot(aes(fill = celltype)) + 
    facet_grid(rows = vars(algorithm)) +
    ggtitle("fitted error parameters") +
    xlab("cell type") +
    ylab("SD(error)")

  # create scatterplots of errors vs offsets
  error_scatter <- ggplot(results.df, aes(x = offset, y = error_sd)) +
    geom_point(aes(col = celltype)) +
    facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
    xlab("Offset") + ylab("Error")

  # create plots of errors + offsets
  plot.df <- data.frame()
  for (algorithm in unique(results.df$algorithm)) {
    for (ct in unique(results.df$celltype)) {
      temp_df <- results.df[
        which(results.df$algorithm == algorithm & results.df$celltype == ct), 
      ]
      offset <- mean(unique(temp_df$offset)) / max(results.df$offset)
      error <- mean(unique(temp_df$error_sd)) / max(results.df$error_sd)
      plot.df <- rbind(
        plot.df,
        data.frame(
          algorithm = algorithm, 
          celltype = ct, 
          combined_score = offset + error
        )
      )
    }
  }
  
  combined_plot <- ggplot(plot.df, aes(x = celltype, y = combined_score)) +
    geom_bar(aes(fill = celltype), stat = "identity") +
    facet_grid(rows = vars(algorithm)) +
    xlab("cell type") +
    ylab("offset + error (scaled)")

  return(list(
    fits = fit_plots,
    offset_plot = offset_plots,
    error_plot = error_plots,
    combined_plot = combined_plot
  ))
}
