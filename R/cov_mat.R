#' create a matrix of covariances between cell type quantities in C
#' 
#' @param C non-negative numeric matrix (cell type x bulk) of cellular 
#' composition 
#' @return matrix of covariances (cell type x cell type)

create_cov_mat <- function(C) {
	cov_mat <- matrix(
    0, nrow = nrow(C), ncol = nrow(C), dimnames = list(rownames(C),rownames(C))
  )
  for (ct in rownames(cov_mat)) {
		for (ct2 in colnames(cov_mat)) {
      cov_mat[ct, ct2] <- cov(C[ct,], C[ct2,])
		}
	}
	return(cov_mat)
}