#' rescale expression profiles to sum to a fixed count number
#'
#' @param exprs expression matrix containing profiles in columns, genes in rows;
#' can be sparse
#' @param count numeric, the number the counts in each profile should sum up to:
#' default NULL, count = nrow(exprs)
#' @return rescaled expression matrix (sparsity depending on type of exprs)

scale_to_count <- function(exprs, count = NULL) {
  # parameter check
  if (!(is.matrix(exprs) || class(exprs) == "dgCMatrix")) {
    stop("exprs must be a matrix")
  }
  # scale total sum either to given value (count)
  # or to the number of rows
  if (is.null(count)) {
    count <- nrow(exprs)
  }
  if (is.matrix(exprs)) {
    csums <- colSums(exprs)
    return(sweep(exprs, 2, csums / count, "/"))
  }else{
    # for sparse matrices, the following is far more efficient
    exprs@x <- unlist(
      lapply(1:(length(exprs@p) - 1),
             function(i, mat = exprs) {
               (mat@x[(mat@p[i] + 1):mat@p[i + 1]] /
               sum(mat@x[(mat@p[i] + 1):mat@p[i + 1]])) * count
             }
      )
    )
    return(exprs)
  }
}
