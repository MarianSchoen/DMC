#' extract matrices of averaged offsets and errors from
#' score simulation results
#' 
#' @param score.results list of data frames as read from h5 result file; 
#' each data frame is one output of `score_algorithms`
#' @return list with two entries:\cr
#' 1) offsets: matrix of offsets (algorithm x cell type)\cr
#' 2) errors: matrix of errors (algorithm x cell type)

create_score_tables <- function (score.results) {
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

  algorithms <- unique(results.df$algorithm)
  celltypes <- unique(results.df$celltype)
  offsets <- matrix(
    0, 
    ncol = length(celltypes),
    nrow = length(algorithms),
    dimnames = list(algorithms, celltypes)
  )
  errors <- matrix(
    0, 
    ncol = length(celltypes),
    nrow = length(algorithms),
    dimnames = list(algorithms, celltypes)
  )
  for (algo in algorithms) {
    for (ct in celltypes) {
      temp_df <- results.df[
        which(results.df$algorithm == algo & results.df$celltype == ct),
      ]
      offsets[algo, ct] <- mean(unique(temp_df$offset))
      errors[algo, ct] <- mean(unique(temp_df$error_sd))
    }
  }
  return(list(offsets = offsets, errors = errors))
}
