#' this functions performs bootstrapping on the real bulks 
#' to estimate the error of the deconvolution results
#' 
#' @param training.exprs matrix containing single-cell expression profiles (training set, one cell per column)
#' @param training.pheno data frame containing phenotype data of the single-cell training set. Has to contain column "cell_type"
#' @param input.algorithms list containing a list for each algorithm. 
#' Each sublist contains 1) name  and 2) function
#' @param verbose logical, default FALSE
#' @param split.data logical, if TRUE (default) then 10% of the training data will be used for reference profile creation and
#' the rest for feature selection/optimization
#' @param exclude.from.bulks character vector containing cell types to be excluded from the bulks (if they are not supplied).
#' If not specified, all will be used.
#' @param exclude.from.signature character vector containing cell types to be excluded from the signature matrix.
#' If not specified, all will be used.
#' @param optimize logical, should the signature matrix be optimized by condition number? If FALSE, max.genes genes will be used
#' @param max.genes maximum number of genes that will be included in the signature for each celltype
#' @param bulks list with two entries:
#' 1) bulks - matrix containing expression data of the bulks (one bulk per column)
#' 2) props - matrix containing the true fractions of cell types within the bulks (cell type x bulk)
#' @return list containing a vector of scores from deconvolution of the bootstrap-samples for each algorithm

bootstrap_bulks <- function(training.exprs, 
                            training.pheno, 
                            algorithms, 
                            verbose = F, 
                            split.data = TRUE, 
                            exclude.from.bulks = NULL, 
                            exclude.from.signature = NULL, 
                            optimize = TRUE, 
                            max.genes = NULL,
                            bulks,
                            props) {
  n.bulks <- ncol(bulks$bulks)
  estimates <- props$est
  real.props <- props$real
  
  cts <- intersect(rownames(estimates[[1]]), rownames(real.props))

  bootstrap.df <- c()
  for(i in 1:20){
    bootstrap.samples <- sample(1:n.bulks, n.bulks, replace = T)
    bootstrap.estimates <- list()
    for(a in names(estimates)){
      bootstrap.estimates[[a]] <- estimates[[a]][, bootstrap.samples]
    }
    bootstrap.real <- real.props[, bootstrap.samples]
    
    for(a in names(estimates)){
      cors <- c()
      for(t in cts){
        temp.cor <- cor(bootstrap.estimates[[a]][t,], bootstrap.real[t,])
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