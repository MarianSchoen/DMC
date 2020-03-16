# written by Tim Mirus

#' evaluate deconvolution performance using different metrics
#' 
#' @param real matrix containing real quantities of cell types (cell type x proportion)
#' @param estimate matrix containing estimations of cell type quantities (cell type x proportion)
#' @return list containing a score for each metric: correlation (cor), mad (mean absolute difference)
#' and rmsd (root mean squared distance)

evaluate_deconvolution <- function(real, estimate) {
  real <- as.vector(real)
  estimate <- as.vector(estimate)

  if (length(real) != length(estimate)) {
    stop("Lengths of quantity vectors do not match")
  }
  correlation <- cor(real, estimate)
  mad <- mean(abs(real - estimate))
  rmsd <- sqrt(mean(abs(real - estimate)^2))
  return(list(cor = correlation, mad = mad, rmsd = rmsd))
}
