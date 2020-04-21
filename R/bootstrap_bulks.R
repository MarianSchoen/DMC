bootstrap_bulks <- function(training.exprs, 
                            training.pheno, 
                            algorithms, 
                            verbose = F, 
                            split.data = TRUE, 
                            exclude.from.bulks = NULL, 
                            exclude.from.signature = NULL, 
                            optimize = TRUE, 
                            max.genes = NULL, 
                            n.bulks = 0, 
                            bulks) {
  n.bulks <- ncol(bulks$bulks)
  bootstrap.results <- data.frame()
  models = NULL
  for(i in 1:50){
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
  colnames(bootstrap.results) <- colnames(deconv.results)
  bootstrap.lst <- list()
  for(a in unique(bootstrap.results$algorithm)){
    bootstrap.lst[[a]] <- bootstrap.results[which(bootstrap.results$algorithm == a), "score"]
  }
  return(bootstrap.lst)
}