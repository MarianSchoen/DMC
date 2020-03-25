#' rescale expression profiles to sum to a fixed count number
#'
#' @param exprs expression matrix containing profiles in columns
#' @param count numeric, the number the counts in each profile should sum up to:
#' default NULL, count = nrow(exprs)
#' @return rescaled expression matrix

scale_to_count <- function(exprs, count = NULL) {
  if (is.null(count)) {
    exprs <- apply(exprs, 2, function(x) {
      (x / sum(x)) * length(x)
    })
  } else {
    exprs <- apply(exprs, 2, function(x) {
      (x / sum(x)) * count
    })
  }
  # remove NAs; they occur if a gene is not expressed in any sample
  if (any(is.na(exprs))) {
    exprs[is.na(exprs)] <- 0
  }
  return(exprs)
}
