#' this functions performs bootstrapping on the real bulk deconvolution results
#' to determine the stability of the deconvolution results
#'
#' @param props list with two entries:
#' 1) est - matrix containing the estimated fractions of cell types
#' within the bulks (cell type x bulk)
#' 2) real - matrix containing the true fractions of cell types
#' within the bulks (cell type x bulk)
#' @param metric evaluation metric; either "cor" (default) or
#' a function conforming to the required format. See benchmark() for details.
#' @return matrix containing all bootstrap runs with columns
#' 'algorithm', 'cell_type' and 'score'

bootstrap_bulks <- function(props, metric = "cor") {
  # parameter check
  if (!is.list(props) || length(props) != 2 ||
      !all(c("est", "real") %in% names(props))) {
        stop("Invalid estimated proportions ('props')")
  }
  if (is.character(metric)) {
    if (metric != "cor") {
      stop("metric must be either \"cor\" or a function")
    }else{
      metric <- cor
    }
  }else{
    if (!is.function(metric)) {
      stop("Function corresponding to 'metric' could not be found.")
    }
  }

  estimates <- props$est
  real.props <- props$real
  n.bulks <- ncol(estimates[[1]])

  cts <- rownames(real.props)
  for (i in seq_len(length(estimates))) {
    cts <- intersect(cts, rownames(estimates[[i]]))
  }
  bootstrap.mat <- c()
  for (i in seq_len(1000)) {
    # draw n.bulks bulks randomly with replacement
    bootstrap.samples <- sample(1:n.bulks, n.bulks, replace = T)

    # create new estimate and real proportion matrix containing these bulks
    bootstrap.estimates <- list()
    for (a in names(estimates)) {
      bootstrap.estimates[[a]] <- estimates[[a]][, bootstrap.samples]
    }
    bootstrap.real <- real.props[, bootstrap.samples, drop = FALSE]

    # calculate for each algorithm for each cell type the correlation
    # between real and estimated proportions
    for (a in names(estimates)) {
      scores <- c()
      for (t in cts) {
        temp.score <- metric(bootstrap.estimates[[a]][t, ], bootstrap.real[t, ])
        # NAs and negative correlations are set to 0
        if (is.na(temp.score) | temp.score < 0) {
          temp.score <- 0
        }
        bootstrap.mat <- rbind(bootstrap.mat, c(a, t, temp.score))
        scores <- c(scores, temp.score)
      }
      score <- mean(scores)
      bootstrap.mat <- rbind(bootstrap.mat, c(a, "overall", score))
    }
  }
  colnames(bootstrap.mat) <- c("algorithm", "cell_type", "score")
  return(bootstrap.mat)
}
