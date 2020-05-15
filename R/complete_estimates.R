#' this function adds rows filled with for all missing cell types 
#' to the estimates returned by an algorithm
#' 
#' @param est.props matrix containing estimated proportions for each cell type (cell type x bulk)
#' @param include.in.x vector containing the cell types that should be present in est.props
#' @return matrix (cell type x bulks) containing the information of est.props

complete_estimates <- function(est.props, include.in.x){
  # check parameters
  if(!is.matrix(est.props)){
    stop("est.props must be a matrix")
  }
  if(!is.character(include.in.x)){
    stop("include.in.x must be a character vector")
  }
  # if any celltype in include.in.x is not in est.props,
  # add a row of zeros with its name as rowname
  for(ct in include.in.x){
    if(!ct %in% rownames(est.props)){
      rnames <- rownames(est.props)
      est.props <- rbind(est.props, rep(NA, ncol(est.props)))
      rownames(est.props) <- c(rnames, ct)
    }
  }
  return(est.props)
}