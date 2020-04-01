# this function adds empty rows for all missing cell types to the estimates returned by an algorithm
complete_estimates <- function(est.props, include.in.x){
  if(!is.matrix(est.props)){
    stop("est.props must be a matrix")
  }
  if(!is.character(include.in.x)){
    stop("include.in.x must be a character vector")
  }
  for(ct in include.in.x){
    if(!ct %in% rownames(est.props)){
      rnames <- rownames(est.props)
      est.props <- rbind(est.props, rep(NA, ncol(est.props)))
      rownames(est.props) <- c(rnames, ct)
    }
  }
  return(est.props)
}