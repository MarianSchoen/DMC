#' cost function to minimize for fitting curves to correlations for different
#' datasets
#' 
#' @param p named vector of parameters that are optimized
#' @param temp_df data frame containing columns 'celltype', 'error', 
#' 'correlation'
#' @param ct character, cell type for which curve should be fit
#' @return residual (loss to be minimized)

cost_fun <- function(
  p,
  temp_df, ct
) {
  sub_df <- temp_df[temp_df$celltype == ct, ]
	e <- p["E"]
	offset <- p["offset"]
  # weights for fit
  if ("error" %in% colnames(sub_df)) {
    w <- sub_df$error
    if(any(w == 0) || any(is.na(w)) || any(is.nan(w))) {
      w <- rep(1, nrow(sub_df))
    }
  }else{
    w <- rep(1, nrow(sub_df))
  }

  prediction <- predict_value(sub_df, offset, e) 
  residual  <- sum(((sub_df$correlation - prediction) / w)^2)
  return(residual)
}
