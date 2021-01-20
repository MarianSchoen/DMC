deconv_for_score <- function(
  training.exprs, training.pheno, 
  test.exprs, test.pheno,
  algorithm.list,
  n.repeats = 3, n.bulks = 100, 
  celltype.column = "celltype", patient.column = "patient", 
  celltypes = NULL, predict.types = NULL, 
  cibersort.path=NULL,
  temp.dir = "./temp", file.suffix = "",
  var.step.size = 5, verbose = TRUE
){
  library(DAB)
  library(Matrix)
  if(!is.null(cibersort.path)){
    source(cibersort.path)
  }
  if(!dir.exists(temp.dir)){
    dir.create(temp.dir)
  }
  
  if(is.null(celltypes)){
    celltypes <- unique(training.pheno[[celltype.column]])
  }
  if(is.null(predict.types)){
    predict.types <- celltypes
  }
  if(!all(celltypes %in% predict.types)){
    excluded.types <- celltypes[which(!(celltypes %in% predict.types))]
  }else{
    excluded.types <- NULL
  }
  
  gc()
  
  ########
  # create bulks and deconvolute them
  results <- list()
  bulk.dir <- paste0(temp.dir, "/bulks")
  if(!dir.exists(bulk.dir)){
    dir.create(bulk.dir)
  }
  for(n in 1:n.repeats){
    # create bulks
    if(!file.exists(paste0(bulk.dir,"/bulks", file.suffix, "_", n, ".RDS"))){
      if(verbose) cat("Create test bulks for deconvolution (Set", n,")...\t", 
                      as.character(Sys.time()), "\n", sep = ""
      )
      bulk.list <- list()
      for(i in c(2^(0:6), 100)){
        cat(i, "\t", sep = "")
        bulk.list[[as.character(i)]] <- create_bulks(
          exprs = test.exprs,
          pheno = test.pheno,
          cell.type.column = celltype.column,
          n.bulks = n.bulks,
          n.profiles.per.bulk = 500,
          sum.to.count = TRUE,
          frac = i
        )
      }
      if(verbose) {
        cat("\n")
        cat("Created bulks successfully. Created list of size ", 
            as.numeric(object.size(bulk.list)) / 1000000,
            " MB\n",
            sep = ""
        )
      }
      saveRDS(bulk.list, paste0(bulk.dir,"/bulks", file.suffix, "_", n, ".RDS"))
    }else{
      bulk.list <- readRDS(paste0(bulk.dir,"/bulks", file.suffix, "_", n, ".RDS"))
    }
    
    # deconvolute using different algorithms
    result.list <- list()
    if(verbose) cat("Deconvolution...\t", as.character(Sys.time()), "\n")
    for(i in names(bulk.list)){
      if(verbose) cat(i, "\t")
      # choose bulks
      bulk.data <- bulk.list[[i]]
      
      result.list[[as.character(i)]] <- list()
      deconv.result <- DAB:::deconvolute(
        training.expr = training.exprs,
        training.pheno = training.pheno,
        test.expr = NULL,
        test.pheno = NULL,
        algorithms = algorithm.list,
        verbose = FALSE,
        exclude.from.signature = excluded.types,
        bulks = bulk.data,
        n.repeats = 1,
        cell.type.column = celltype.column,
        patient.column = patient.column	
      )$results.list
      for(name in names(deconv.result[[1]])){
        result.list[[as.character(i)]][name] <- deconv.result[[1]][[name]]$est.props
      }
      if(verbose) cat("\n")
    }
    results[[n]] <- result.list
  }
  results.dir <- paste0(temp.dir, "/results")
  if(!dir.exists(results.dir)){
    dir.create(results.dir)
  }
  saveRDS(results, paste0(results.dir, "/deconv_results.RDS"))
  return(results)
}

process_score_deconv_results <- function(results.all, celltypes, temp.dir = "./temp", verbose = TRUE)
{
  library(ggplot2)
  library(DAB)
  library(Matrix)
  
  if(verbose) cat("Calculating correlations...\t", as.character(Sys.time()), "\n")
  
  cor.df <- c()
  for(n in 1:length(results.all)){
    results <- results.all[[n]]
    for(i in names(results)){
      cors <- list()
      for(name in names(results[[i]])){
        if(name == "real.props") next
        if(is.matrix(results[[i]][[name]])){
          cors[[name]] <- sapply(celltypes, 
                                 function(x, props1 = results[[i]][[name]], props2 = results[[i]]$real.prop){
                                   cor(props1[x,], props2[x,])
                                 })
        }else{
          cors[[name]] <- rep(0, length(celltypes))
        }
        names(cors[[name]]) <- celltypes
      }
      
      for(ct in celltypes){
        ct.sd <- sd(results[[i]]$real.props[ct,])
        for(name in names(cors)){
          cor.df <- rbind(cor.df, c(ct, cors[[name]][ct], name, i, ct.sd, n))
        }
      }
    }
  }
  cor.df <- as.data.frame(cor.df)
  colnames(cor.df) <- c("celltype", "score", "algorithm", "index", "ct_sd", "repetition")
  cor.df$score <- as.numeric(as.character(cor.df$score))
  cor.df$index <- as.numeric(as.character(cor.df$index))
  cor.df$ct_sd <- as.numeric(as.character(cor.df$ct_sd))
  cor.df$repetition <- as.numeric(as.character(cor.df$repetition))
  
  if(any(is.na(cor.df$score))){
    cor.df[which(is.na(cor.df$score)), "score"] <- 0
  }
  
  results.dir <- paste0(temp.dir, "/results")
  if(!dir.exists(results.dir)){
    dir.create(results.dir)
  }
  saveRDS(cor.df, paste0(results.dir, "/correlations.RDS"))
  return(cor.df)
}

fit_model <- function(correlations, temp.dir = "./temp", verbose = TRUE){
  library(DAB)
  library(ggplot2)
  
  results.dir <- paste0(temp.dir, "/results")
  if(!dir.exists(results.dir)){
    dir.create(results.dir)
  }
  
  plot.dir <- paste0(temp.dir, "/plots")
  if(!dir.exists(plot.dir)){
    dir.create(plot.dir)
  }
  
  # reduce to data frame with mean and sd across repetitions
  reduced.correlations <- data.frame()
  
  if(verbose) cat("Calculate mean and sd of correlation for each algorithm, cell type and variance...\n")
  for(ct in unique(correlations$celltype)){
    for(algo in unique(correlations$algorithm)){
      for(i in unique(correlations$index)){
        selection <- intersect(
          which(correlations$celltype == ct),
          intersect(
            which(correlations$algorithm == algo), 
            which(correlations$index == i)
          )
        )
        m <- mean(correlations[selection,"score"])
        e <- sd(correlations[selection, "score"])
        ct_sd <- mean(correlations[selection, "ct_sd"])
        
        reduced.correlations <- rbind(
          reduced.correlations,
          data.frame(algorithm = algo, index = i, celltype = ct, score = m, sd = e, ct_sd = ct_sd)
        )
      }
    }
  }
  correlations <- reduced.correlations
  rm(reduced.correlations)
  saveRDS(correlations, paste0(results.dir, "/correlations_sd.RDS"))
  
  cor.mat <- matrix(NA, ncol = length(unique(correlations$algorithm)), nrow = length(unique(correlations$celltype)))
  colnames(cor.mat) <- unique(correlations$algorithm)
  rownames(cor.mat) <- unique(correlations$celltype)
  
  overall.fit.df <- data.frame()
  pdf(paste0(plot.dir, "/model_fits.pdf"), width = 8 * length(unique(correlations$algorithm)), height = 8)
  fits <- list()
  for(ct in unique(correlations$celltype)){
    if(verbose) cat(ct, "\n")
    ct.fit.df <- data.frame() # contains results
    algo.res <- c() # algorithm score as string
    algo.names <- c() # algorithm names for algo.res
    fits[[ct]] <- list()
    for(algo in unique(correlations$algorithm)){
      if(verbose) cat(algo, "...\t")
      temp.cors <- correlations[which(correlations$algorithm == algo),]
      temp.cors <- temp.cors[which(temp.cors$celltype == ct),]
      
      # fit a model to the scores
      fitted.model <- try(nls(
        score ~ 1 / sqrt(1 + (x/index)^2),
        data = temp.cors, 
        start = list(x = 1),
        weights = temp.cors$sd^(-2),
        control = nls.control(maxiter = 1000000, tol = 1e-5)
      ))
      if(class(fitted.model) == "try-error"){
        next
      }
      fits[[ct]][[algo]] <- fitted.model
      algo.names <- c(algo.names, algo)
      
      fitted.x <- abs(as.numeric(fitted.model$m$getPars()["x"]))
      
      # create data frame containing values predicted from fit
      fit.df <- data.frame(
        correlation = temp.cors$score, 
        sd = temp.cors$index, 
        prediction = predict(fitted.model, newdata = temp.cors$index), 
        algorithm = algo,
        sd_error = fitted.x,
        score_sd = temp.cors$sd
      )
      
      # one data frame contains all algorithm results for this cell type
      ct.fit.df <- rbind(ct.fit.df, fit.df)
      cor.mat[ct, algo] <- temp.cors[which(temp.cors$index == 8), "score"]
      
      algo.res <- c(algo.res, paste0("SD(error): ", round(fitted.x,4)))
    }
    names(algo.res) <- algo.names
    ct.fit.df$algorithm <- factor(ct.fit.df$algorithm, levels = names(algo.res))
    
    # combine data frames from different cell types
    overall.fit.df <- rbind(overall.fit.df, ct.fit.df)
    
    # plot
    labs = paste(names(algo.res), algo.res, sep = "\n")
    names(labs) <- names(algo.res)
    p <- ggplot(ct.fit.df, aes(x = sd, y = correlation)) + 
      geom_pointrange(aes(ymin = correlation - score_sd, ymax = correlation + score_sd), col = "black") + 
      geom_line(aes(x = sd, y = prediction), col = "red") + 
      ggtitle("calculated correlation and fitted model", subtitle = paste0(ct)) + 
      xlab("bulk data SD") + ylab("Correlation") + ylim(0,1)+
      facet_grid(cols = vars(algorithm), labeller = labeller(algorithm = labs))
    plot(p)
    
    if(verbose) cat("\n")
  }
  dev.off()
  
  # produce a number for each algorithm and plot overall correlations
  algo.scores <- sapply(unique(correlations$algorithm),
                        function(x, df = overall.fit.df){
                          temp.df <- df[which(df$algorithm == x),]
                          return(mean(unique(temp.df$sd_error)))
                        }
  )
  if(verbose) {
    cat("Scores:\n")
    print(algo.scores)
  }
  
  # for each algorithm average over all celltypes
  overall.df <- data.frame()
  for(algo in unique(correlations$algorithm)){
    temp.df <- correlations[which(correlations$algorithm == algo),]
    for(i in unique(correlations$index)){
      overall.df <- rbind(overall.df,
                          data.frame(algorithm = algo, index = i, celltype = "overall", score = mean(temp.df[which(temp.df$index == i),"score"]), sd = mean(temp.df[which(temp.df$index == i),"ct_sd"]))
      )
    }
  }
  
  # plot final algorithm scores with averaged curve
  pdf(paste0(plot.dir,"/final_scores.pdf"), width = 16, height = 9)
  for(algo in unique(overall.df$algorithm)){
    if(verbose) cat(algo, "\n")
    p1 <- ggplot(overall.df[which(overall.df$algorithm == algo),], aes(x = index, y = score)) + 
      ggtitle(paste0("Overall Correlation for ", algo), subtitle = paste0("Score [mean SD(Error)]: ", round(algo.scores[algo],2))) + 
      xlab("sd of cell types in bulks") + ylab("correlation") + geom_line() + ylim(0,1)
    plot(p1)
  }
  dev.off()
  
  # plot algorithm score vs correlation
  cors <- apply(cor.mat, 2, mean)
  cors <- cors[names(algo.scores)]
  
  temp.df <- data.frame(algorithm = names(algo.scores), score = algo.scores, correlation = cors)
  p2 <- ggplot(temp.df, aes(x = correlation, y = score)) + geom_line() + geom_point(aes(col = algorithm)) + ggtitle("Score vs. Correlation for each algorithm") + xlim(0,1)
  pdf(paste0(plot.dir, "/score_vs_cor.pdf"), width = 16, height = 9)
  plot(p2)
  dev.off()
  
  return(list(Scores = algo.scores, fits = fits))
}

score_algorithms <- function(
  training.exprs, training.pheno, 
  test.exprs, test.pheno,
  algorithm.list,
  n.repeats = 3, n.bulks = 100, 
  celltype.column = "celltype", patient.column = "patient", 
  celltypes = NULL, predict.types = NULL, 
  cibersort.path=NULL,
  temp.dir = "./temp", file.suffix = "",
  var.step.size = 5, verbose = TRUE
){
  
  deconv.results <- deconv_for_score(
    training.exprs, training.pheno, 
    test.exprs, test.pheno,
    algorithm.list,
    n.repeats = n.repeats, n.bulks = n.bulks, 
    celltype.column = celltype.column, patient.column = patient.column, 
    celltypes = celltypes, predict.types = predict.types, 
    cibersort.path=cibersort.path,
    temp.dir = temp.dir, file.suffix = file.suffix,
    var.step.size = var.step.size, verbose = verbose)
  
  processed_results <- process_score_deconv_results(deconv.results, celltypes, temp.dir)
  
  results <- fit_model(processed_results, temp.dir, verbose)
  return(results$Scores)
}