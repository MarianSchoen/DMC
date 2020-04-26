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
                            bulks) {
  n.bulks <- ncol(bulks$bulks)
  bootstrap.results <- data.frame()
  models = NULL
  for(i in 1:50){
    if(verbose) cat(i, " ", sep = "")
    bootstrap.samples <- sample(1:n.bulks, n.bulks, replace = T)
    bootstrap.bulks <- list(bulks = bulks$bulks[, bootstrap.samples], props = bulks$props[,bootstrap.samples])
    colnames(bootstrap.bulks$bulks) <- 1:n.bulks
    colnames(bootstrap.bulks$props) <- 1:n.bulks
    deconv.results <- deconvolute(
      training.exprs, 
      training.pheno, 
      NULL, NULL, 
      algorithms, 
      verbose, TRUE, NULL, 
      exclude.from.signature, 
      TRUE, NULL, 0,
      bootstrap.bulks,
      n.repeats = 1,
      algorithm.models = models
    )
    #print(str(deconv.results))
    if(i == 1){
      models <- list()
      for(a in algorithms){
        ref.profiles <- deconv.results$results.list[["1"]][[a$name]]$ref.profiles
        g <- deconv.results$results.list[["1"]][[a$name]]$g
        models[[a$name]] <- list(ref.profiles = ref.profiles, g = g)
      }
    }
    deconv.results <- prepare_data(results.all = deconv.results, metric = "cor")
    bootstrap.results <- rbind(bootstrap.results, deconv.results[which(deconv.results$cell_type == "overall"),])
  }
  cat("\n")
  colnames(bootstrap.results) <- colnames(deconv.results)
  bootstrap.lst <- list()
  for(a in unique(bootstrap.results$algorithm)){
    bootstrap.lst[[a]] <- bootstrap.results[which(bootstrap.results$algorithm == a), "score"]
  }
  return(bootstrap.lst)
}