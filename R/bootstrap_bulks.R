#' this functions performs bootstrapping on the real bulks 
#' to estimate the error of the deconvolution results
#' 
#' @param props list with two entries:  
#' 
#' 1) est - matrix containing the estimated fractions of cell types within the bulks (cell type x bulk)  
#' 
#' 2) real - matrix containing the true fractions of cell types within the bulks (cell type x bulk)
#' @return list containing a vector of scores from deconvolution of the bootstrap-samples for each algorithm

bootstrap_bulks <- function(props) {

  # parameter check
  if(!is.list(props) || length(props) != 2 || !all(c("est", "real") %in% names(props))){
    stop("Invalid estimated proportions ('props')")
  }

  estimates <- props$est
  real.props <- props$real
  n.bulks <- ncol(estimates)
  cts <- intersect(rownames(estimates[[1]]), rownames(real.props))

  bootstrap.df <- c()
  for(i in 1:1000){
    # draw n.bulks bulks randomly with replacement
    bootstrap.samples <- sample(1:n.bulks, n.bulks, replace = T)

    # create new estimate and real proportion matrix containing these bulks
    bootstrap.estimates <- list()
    for(a in names(estimates)){
      bootstrap.estimates[[a]] <- estimates[[a]][, bootstrap.samples]
    }
    bootstrap.real <- real.props[, bootstrap.samples]
    
    # calculate for each algorithm for each cell type the correlation 
    # between real and estimated proportions
    for(a in names(estimates)){
      cors <- c()
      for(t in cts){
        temp.cor <- cor(bootstrap.estimates[[a]][t,], bootstrap.real[t,])
        # NAs and negative correlations are set to 0
        if(is.na(temp.cor) | temp.cor < 0){
          temp.cor <- 0
        }
        bootstrap.df <- rbind(bootstrap.df, c(a, t, temp.cor))
        cors <- c(cors, temp.cor)
      }
      score <- mean(cors)
      bootstrap.df <- rbind(bootstrap.df, c(a, 'overall', score))
    }
  }
  colnames(bootstrap.df) <- c("algorithm", "cell_type", "score")
  return(bootstrap.df)
}