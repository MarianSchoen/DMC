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
  # reduce repetitions by averaging
  reduced_df <- data.frame()
  for (a in results.df$algorithm) {
    for (cs in results.df$cluster_size) {
      idx <- which(
        results.df$algorithm == a &
        results.df$cluster_size == cs
      )
      idx_fine <- intersect(idx, which(!results.df$coarse))
      idx_coarse <- intersect(idx, which(results.df$coarse))
      
      avg_score_fine <- mean(results.df[idx_fine, "score"])
      avg_score_coarse <- mean(results.df[idx_coarse, "score"])
      
      sd_score_fine <- sd(results.df[idx_fine, "score"])
      sd_score_coarse <- sd(results.df[idx_coarse, "score"])
      reduced_df <- rbind(
        reduced_df, 
        data.frame(
          algorithm = a,
          cluster_size = cs,
          coarse = FALSE,
          score = avg_score_fine,
          error = sd_score_fine
        ),
        data.frame(
          algorithm = a,
          cluster_size = cs,
          coarse = TRUE,
          score = avg_score_coarse,
          error = sd_score_coarse
        )
      )
    }
  }
  
  # plotting
  overall.plot <- ggplot(
    reduced_df[which(!reduced_df$coarse), ],
    aes(x = cluster_size, y = score, col = algorithm)
  ) +
  geom_line() +
  geom_errorbar(
    aes(ymin = score - error, ymax = score + error, col = algorithm)
  ) +
  scale_x_reverse()

  algorithm.plot <- ggplot(
    reduced_df,
    aes(x = cluster_size, y = score, col = algorithm)
  ) +
  geom_line(aes(linetype = coarse)) +
  geom_errorbar(
    aes(ymin = score - error, ymax = score + error, col = algorithm)
  ) +
  scale_x_reverse() +
  facet_grid(rows = vars(algorithm))

  reduced_df$coarse <- ifelse(as.logical(reduced_df$coarse), "coarse", "fine")

  coarse.plot <- ggplot(
    reduced_df,
    aes(x = cluster_size, y = score, col = algorithm)
  ) +
  scale_x_reverse() +
  geom_line() +
  geom_errorbar(
    aes(ymin = score - error, ymax = score + error, col = algorithm)
  ) +
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
