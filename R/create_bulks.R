#' simulate bulk data by summing over single-cell data
#'
#' @param exprs non-negative numeric, scRNA-seq profiles as columns, 
#' features as rows
#' @param pheno data.frame, phenotype data containing
#' cell type labels for the expression matrix,
#' must contain 'cell.type.column',
#' ordering of cells must be the same as in exprs
#' @param n.bulks integer, the number of bulks to be created, defaults to 500
#' @param include.in.bulks vector of strings, cell types to be used for bulk 
#' simulation; if not supplied, all will be used
#' @param n.profiles.per.bulk positive numeric, number of samples to be randomly
#' drawn for each bulk; default 1000
#' @param sum.to.count boolean, should all bulks be normalized
#' to a fixed total count number? default TRUE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? 
#' @param frac numeric >= 1; determines the variance of cell type proportions within the simulated bulks;
#' higher d means more extreme distributions
#' @param switch.prob fraction of bulks sampled with inverse cell type quantity
#' @return list with 
#'    - "bulks": matrix containing bulk expression profiles (features x bulks)
#'    - "props" matrix containing quantities (cell type x bulks)
#' @example create_bulks(training.exprs, training.pheno, n.bulks = 1000)
#' @export

create_bulks <- function(
  exprs, 
  pheno, 
  cell.type.column = "cell_type",
  n.bulks = 500, 
  include.in.bulks = NULL, 
  n.profiles.per.bulk = 1000, 
  sum.to.count = TRUE,
  frac = 5,
  switch.prob = 0
) {
  # parameter checks
  if (nrow(pheno) != ncol(exprs)) {
    stop("Number of columns in exprs and rows in pheno do not match")
  }
  if(!is.numeric(n.bulks)){
    stop("n.bulks must be numeric")
  }
  if(!is.numeric(n.profiles.per.bulk)){
    stop("n.profiles.per.bulk must be numeric > 0")
  }else{
    if(n.profiles.per.bulk <= 0){
      stop("n.profiles.per.bulk must be numeric > 0")
    }
  }
  rownames(pheno) <- colnames(exprs)
  
  # keep only specified cell types
  if (is.null(include.in.bulks)) {
    include.in.bulks <- unique(pheno[, cell.type.column])
  }
  to.remove <- which(pheno[, cell.type.column] == "unassigned")
  if(length(to.remove) > 0){
    pheno <- pheno[-to.remove,]
    exprs <- exprs[,-to.remove]
  }
  
  if (length(which(pheno[, cell.type.column] %in% include.in.bulks)) > 0) {
    pheno <- pheno[which(pheno[, cell.type.column] %in% include.in.bulks), , drop = FALSE]
  }
  exprs <- exprs[, rownames(pheno), drop = FALSE]
  
  # create a matrix to contain the bulk expression profiles
  bulk.exprs <- matrix(
    0,
    nrow = nrow(exprs),
    ncol = n.bulks
  )
  rownames(bulk.exprs) <- rownames(exprs)
  colnames(bulk.exprs) <- as.character(1:n.bulks)
  
  # create a matrix to contain true proportions for each simulated bulk
  props <- matrix(0, nrow = length(unique(pheno[, cell.type.column])), ncol = n.bulks)
  rownames(props) <- unique(pheno[, cell.type.column])
  colnames(props) <- colnames(bulk.exprs)
  
  types <- unique(pheno[, cell.type.column])
  types <- sample(types, length(types), replace = FALSE)
  
  # sample random weights for each type	
  ct.fracs <- sort(table(pheno[, cell.type.column]) / nrow(pheno) * 100, decreasing = TRUE)

  if(switch.prob > 0){
  	ct.fracs.2 <- 100 - ct.fracs
  	names(ct.fracs.2) <- types
  }else{
	  ct.fracs.2 <- ct.fracs
  }

  for(i in seq_len(n.bulks)){
    bulk.samples <- c()
    weights <- c()
    if (runif(1, 0, 1) < switch.prob){
	    switch <- TRUE
    }else{
	    switch <- FALSE
    }
    for(ct in types){
      if(switch){
	      ct.frac <- ct.fracs[ct]
      }else{
	      ct.frac <- ct.fracs.2[ct]
      }

      if(frac / 2 > ct.frac){
        frac.min <- 1
        frac.max <- frac
      }else if(frac / 2 > 100 - ct.frac){
        frac.max <- 100
        frac.min <- 100 - frac
      }else{
        frac.min <- ct.frac - frac / 2
        frac.max <- ct.frac + frac / 2
      }
      frac.min <- as.integer(frac.min)
      frac.max <- as.integer(frac.max)
      w <-  sample(frac.min:frac.max, 1)
      
      if(w < 1) w <- 1
      if(w > 100) w <- 100
      weights <- c(weights, w)
    }
    names(weights) <- types
    if(any(weights < 0)){
      weights[weights < 0] <- 0
    }
    weights <- weights/sum(weights)

    
    # determine, how many samples of each type should be drawn
    n.samples <- sample(types, size = ceiling(n.profiles.per.bulk), replace = TRUE, prob = weights) 
    type.table <- table(n.samples)
    
    # draw samples for each type, store the proportions and the expression
    for(t in types){
      if(t %in% names(type.table)){
        coarse.samples <- sample(which(pheno[,cell.type.column] == t), size = type.table[t], replace = TRUE)
        bulk.samples <- c(bulk.samples, coarse.samples)
        props[t, i] <- length(coarse.samples)
      }
    }
    props[,i] <- props[,i] / length(bulk.samples)
    
    bulk.exprs[,i] <- Matrix::rowSums(exprs[, bulk.samples, drop = F])
  }
  
  # let no feature occur twice in the bulks
  if(any(duplicated(rownames(bulk.exprs)))){
    warning("Found duplicate features in simulated bulks. Removing...")
    bulk.exprs <- bulk.exprs[-which(duplicated(bulk.exprs)), ]
  }
  
  # sum bulks to fixed total count per profile if sum.to.count is true
  if (sum.to.count) {
    bulk.exprs <- scale_to_count(bulk.exprs)
  }
  
  return(list(bulks = bulk.exprs, props = props))
}