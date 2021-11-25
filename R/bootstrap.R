bootstrap_results <- function(real, estimate, nrep = 1000) {
  if (!is.numeric(real) || !is.numeric(estimate)) {
    warning("Parameters 'real' and 'estimate' need to be numeric vectors")
    return(NULL)
  }
  if (length(real) != length(estimate)) {
    warning("Length of real and estimated quantities does not match")
    return(NULL)
  }
  correlations <- c()
  for (i in 1:nrep) {
    samples <- sample(1:length(estimate), length(estimate), replace = TRUE)
    correlations <- c(correlations, cor(real[samples], estimate[samples]))
  }
  #print(cor(real, estimate))
  #print(mean(correlations))
  return(list(mean_cor = mean(correlations), error = sd(correlations)))
}
