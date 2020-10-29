#' rescale expression profiles to sum to a fixed count number
#'
#' @param exprs expression matrix containing profiles in columns
#' @param count numeric, the number the counts in each profile should sum up to:
#' default NULL, count = nrow(exprs)
#' @return rescaled expression matrix

scale_to_count <- function(exprs, count = NULL) {
  # parameter check
  if(!is.matrix(exprs)){
    stop("exprs must be a matrix")
  }
  # scale total sum either to given value (count)
  # or to the number of rows 
  if (is.null(count)) {
    count <- nrow(exprs)
  }
  csums <- colSums(exprs)
  return(sweep(exprs, 2, csums / count, "/"))

  # remove NAs; they occur if a gene is not expressed in any sample
  #if (any(is.na(exprs))) {
  #  exprs[is.na(exprs)] <- 0
  #}
  #return(exprs)
}
