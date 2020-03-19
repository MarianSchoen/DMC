#' via_dtd_get_random_data 
#' 
#' Generate random cell type like data, and generate bulks from them
#'
#' @param n.features integer, number of features
#' @param n.bulks integer, number of bulks
#' @param n.types.in.data integer, number of simulated cell types in the data
#' @param n.types.to.deconvolute integer, number of simulated cell types
#' @param per.type integer, profiles per cell type in the data
#'
#' @return named list with 4 entries (all numeric matrices)
#' @export
via_dtd_get_random_data <- function(
  n.features = 100, 
  n.bulks = 100, 
  n.types.in.data = 20,
  n.types.to.deconvolute = 10,
  per.type = 30
){
  
  
  # I only check for numeric, as integers are not really integers in R
  if(any(
    !is.numeric(n.features) ||
    !is.numeric(n.bulks) ||
    !is.numeric(n.bulks) ||
    !is.numeric(n.types.in.data) ||
    !is.numeric(n.types.to.deconvolute) ||
    !is.numeric(per.type)
  )){
    stop("In via_dtd_get_random_data: all input parameters must be integers")
  }
  
  
  if(n.types.to.deconvolute > n.types.in.data){
    warning(
    "You can't deconvolute more types than there are in the data 
    => therefore, resetting 'n.types.to.deconvolute'"
    )
    n.types.to.deconvolute <- n.types.in.data
  }
  
  library(DTD)
  set.seed(1)
  
  some.data <- generate_random_data(
    n.types = n.types.in.data
    , n.samples.per.type = per.type
    , n.features = n.features
  )
  
  indi.vec <- gsub(
    pattern = "^.*.T"
    , replacement = "T"
    , x = colnames(some.data)
  )
  names(indi.vec) <- colnames(some.data)
  
  bulk.matrix <- matrix(ncol = 0, nrow = n.features)
  rownames(bulk.matrix) <- rownames(some.data)
  
  c.matrix <- matrix(ncol = 0, nrow = n.types.to.deconvolute)
  rownames(c.matrix) <- paste0("Type", 1:n.types.to.deconvolute)

  for(random.bulk in 1:n.bulks){
    pos <- sample(
      x = 1:ncol(some.data)
      , size = n.features # false friend, but as long as it works ...
    )
    
    mixture <- rowSums(some.data[, pos])
    bulk.matrix <- cbind(bulk.matrix, mixture)
    
    c.vec <- table(indi.vec[pos])
    # check wheter all types are in 'c.vec'
    missing.types <- which(!rownames(c.matrix) %in% names(c.vec))
    if(length(missing.types) > 0){
      missing <- paste0("Type", missing.types)
      missing.vec <- rep(0, length(missing))
      names(missing.vec) <- missing
      c.vec <- c(c.vec, missing.vec)
    }
    c.matrix <- cbind(c.matrix, c.vec[rownames(c.matrix)])
  }
  colnames(c.matrix) <- paste0("mixture", 1:ncol(c.matrix))
  colnames(bulk.matrix) <- paste0("mixture", 1:ncol(bulk.matrix))
  
  ret <- list(
    "sc.data" = some.data, 
    "sc.pheno" = indi.vec, 
    "bulks" = bulk.matrix, 
    "bulks.pheno" = c.matrix
  )
  return(ret)
}