#' create lineplot of algorithm performance across different cell type
#' granularities
#'
#' @param results.df data frame as returned by prepare_data
#' @return list containing\cr
#' 1) lineplot of all algorithms\cr
#' 2) lineplots of individual algorithms\cr
#' 3) lineplot of coarse correlations after fine deconvolution

subtype_plot <- function(
  results.df
) {
  results.df <- results.df[which(results.df$cell_type == "overall"), ]
  overall.plot <- ggplot(
    results.df[which(!as.logical(results.df$coarse)), ],
    aes(x = cluster_size, y = score, col = algorithm)
  ) +
  geom_line() +
  scale_x_reverse()

  algorithm.plot <- ggplot(
    results.df,
    aes(x = cluster_size, y = score, col = algorithm)
  ) +
  geom_line(aes(linetype = coarse)) +
  scale_x_reverse() +
  facet_grid(rows = vars(algorithm))

  results.df$coarse <- ifelse(as.logical(results.df$coarse), "coarse", "fine")

  coarse.plot <- ggplot(
    results.df,
    aes(x = cluster_size, y = score, col = algorithm)
  ) +
  scale_x_reverse() +
  geom_line() +
  facet_grid(
    cols = vars(coarse)
  ) +
  theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22),
      title = element_text(size = 24),
      axis.title.x = element_text(size = 22),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 20),
      strip.text.x = element_text(size = 20)
    ) +
    xlab("average number of profiles per subtype")

  return(list(
    overall = overall.plot,
    algorithm = algorithm.plot,
    coarse = coarse.plot
  ))
}
