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
        switch.prob <- 0.4
       
	if(i < 40){
		switch.prob <- 0.2
	}
	 if(i < 10){
		switch.prob <- 0
	 }
        bulk.list[[as.character(i)]] <- create_bulks(
          exprs = test.exprs,
          pheno = test.pheno,
          cell.type.column = celltype.column,
          n.bulks = n.bulks,
          include.in.bulks = celltypes,
          n.profiles.per.bulk = 500,
          sum.to.count = TRUE,
          frac = i,
          switch.prob = switch.prob
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
      cat("\n\n\n")
      for(name in names(deconv.result[[1]])){
        result.list[[as.character(i)]][[name]] <- deconv.result[[1]][[name]]$est.props
      }
      result.list[[as.character(i)]][["real.props"]] <- bulk.data$props
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
                                 function(x, props1 = results[[i]][[name]], props2 = results[[i]]$real.props){
					 if((x %in% rownames(props1)) && (x %in% rownames(props2))){
                                   return(cor(props1[x,], props2[x,]))
					 }else{
						 return(NA)
					 }
                                 })
        }else{
          cors[[name]] <- rep(0, length(celltypes))
        }
        names(cors[[name]]) <- celltypes
        if(any(is.na(cors[[name]]))){
          cors[[name]][is.na(cors[[name]])] <- 0
        }
        if(any(cors[[name]] < 0)){
          cors[[name]][cors[[name]] < 0] <- 0
        }
      }
      
      
      # calculate bulk entropy here
      bulk.entropies <- apply(
        results[[i]]$real.props, 
        2, 
        DescTools::Entropy
      )
      # normalize by maximum entropy
      final.entropy <- sum(bulk.entropies) / (ncol(results[[i]]$real.props) * log2(results[[i]]$real.props))
      
      for(ct in celltypes){
	      if(ct %in% rownames(results[[i]]$real.props)){
        	ct.sd <- sd(results[[i]]$real.props[ct,])
	      }else{
		      ct.sd <- 0
	      }
        for(name in names(cors)){
          cor.df <- rbind(cor.df, c(ct, cors[[name]][ct], name, i, ct.sd, n, final.entropy))
        }
      }
    }
  }
  cor.df <- as.data.frame(cor.df)
  colnames(cor.df) <- c("celltype", "score", "algorithm", "index", "ct_sd", "repetition", "entropy")
  cor.df$score <- as.numeric(as.character(cor.df$score))
  cor.df$index <- as.numeric(as.character(cor.df$index))
  cor.df$ct_sd <- as.numeric(as.character(cor.df$ct_sd))
  cor.df$entropy <- as.numeric(as.character(cor.df$entropy))
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

fit_model <- function(correlations, temp.dir = "./temp", verbose = TRUE, repetitions = 1){
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
        entropy <- mean(correlations[selection, "entropy"])
        
        reduced.correlations <- rbind(
          reduced.correlations,
          data.frame(algorithm = algo, index = i, celltype = ct, score = m, sd = e, ct_sd = ct_sd, entropy = entropy)
        )
      }
    }
  }
  correlations <- reduced.correlations
  rm(reduced.correlations)
  saveRDS(correlations, paste0(results.dir, "/correlations_sd.RDS"))
  
  # plot entropy 
  pdf(paste0(plot.dir, "/entropy_plots.pdf"), width = 12, height = 8)
      p <- ggplot(correlations, aes(x = entropy, y = score)) +
        geom_point() + geom_line(aes(col = celltype)) + geom_errorbar(aes(ymin = score - sd, ymax = score + sd), alpha = 0.7) +
        facet_grid(rows = vars(celltype), cols = vars(algorithm))
      plot(p)
      p <- ggplot(correlations, aes(x = index, y = entropy)) + geom_point()
      plot(p)
  dev.off()
  #############
  
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
      
      if(any(is.na(temp.cors$sd))){
	      weights <- rep(1, length(temp.cors$sd))
      }else{
	      weights <- temp.cors$sd^(-2)
      }
      fitted.model <- try(nls(
        score ~ (1 / sqrt(1 + (x/ct_sd)^2)) * alpha,
        data = temp.cors, 
        start = list(x = 0.01, alpha = 0.5),
        weights = weights,
        control = nls.control(maxiter = 1000000, tol = 1e-5),
        algorithm = "port",
        lower = c(0, 0),
        upper = c(Inf, 1)
      ))
      if(class(fitted.model) == "try-error"){
        fits[[ct]][[algo]] <- list(NULL)
        fitted.x <- NA
	      fitted.alpha <- NA
        
        # create data frame containing values predicted from fit
        fit.df <- data.frame(
          correlation = temp.cors$score, 
          sd = temp.cors$ct_sd, 
          prediction = 0, 
          algorithm = algo,
          sd_error = fitted.x,
          offset = fitted.alpha,
          score_sd = temp.cors$sd
        )
        fit.df$sd_error <- as.numeric(fit.df$sd_error)
      }else{
        fits[[ct]][[algo]] <- fitted.model
        fitted.x <- abs(as.numeric(fitted.model$m$getPars()["x"]))
        fitted.alpha <- abs(as.numeric(fitted.model$m$getPars()["alpha"]))
        # create data frame containing values predicted from fit
        fit.df <- data.frame(
          correlation = temp.cors$score, 
          sd = temp.cors$ct_sd, 
          prediction = predict(fitted.model, newdata = temp.cors$ct_sd), 
          algorithm = algo,
          sd_error = fitted.x,
          offset = fitted.alpha,
          score_sd = temp.cors$sd
        )
      }
      
      algo.names <- c(algo.names, algo)
      
      # one data frame contains all algorithm results for this cell type
      ct.fit.df <- rbind(ct.fit.df, fit.df)
      cor.mat[ct, algo] <- temp.cors[which(temp.cors$index == 8), "score"]
      
      # error and offset
      algo.res <- c(algo.res, paste0("SD(error): ", round(fitted.x,4), "\nOffset: ", round(fitted.alpha, 3)))
    }
    names(algo.res) <- algo.names
    ct.fit.df$algorithm <- factor(ct.fit.df$algorithm, levels = names(algo.res))
    
    # combine data frames from different cell types
    overall.fit.df <- rbind(overall.fit.df, ct.fit.df)
    
    if(nrow(ct.fit.df) > 0){
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
    }else{
      labs = paste(names(algo.res), algo.res, sep = "\n")
      names(labs) <- names(algo.res)
      p <- ggplot(ct.fit.df, aes(x = sd, y = correlation)) + 
        geom_pointrange(aes(ymin = correlation - score_sd, ymax = correlation + score_sd), col = "black") + 
        geom_line(aes(x = sd, y = prediction), col = "red") + 
        ggtitle("calculated correlation and fitted model", subtitle = paste0(ct)) + 
        xlab("bulk data SD") + ylab("Correlation") + ylim(0,1)+
        facet_grid(cols = vars(algorithm), labeller = labeller(algorithm = labs))
      plot(p)
    }
    if(verbose) cat("\n")
  }
  dev.off()
  if(nrow(overall.fit.df) == 0){
    warning("All fits failed. Returning NULL.")
    return(NULL)
  }
  
  # produce a number for each algorithm and plot overall correlations
  
  # make this two numbers
  algo.scores <- lapply(unique(correlations$algorithm),
                        function(x, df = overall.fit.df){
                          temp.df <- df[which(df$algorithm == x),]
                          if(any(is.na(temp.df$sd_error))){
                            mean.sd <- mean(unique(temp.df$sd_error[!is.na(temp.df$sd_error)]))
                          }else{
                            mean.sd <- mean(unique(temp.df$sd_error))
                          }
                          
                          if(any(is.na(temp.df$offset))){
                            mean.offset <- mean(unique(temp.df$offset[!is.na(temp.df$offset)]))
                          }else{
                            mean.offset <- mean(unique(temp.df$offset))
                          }
                          return(list(error = mean.sd, offset = mean.offset))
                        }
  )
  names(algo.scores) <- unique(correlations$algorithm)

  if(verbose) {
    cat("Scores:\n")
    cat("\t\tError\tOffset\n")
    for(algo in names(algo.scores)){
      cat(algo, ": \t", algo.scores[[algo]]$error,"\t", algo.scores[[algo]]$offset, "\n")
    }
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
    p1 <- ggplot(overall.df[which(overall.df$algorithm == algo),], aes(x = index, y = score)) + 
      ggtitle(paste0("Overall Correlation for ", algo), subtitle = paste0("Score [mean SD(Error)]: ", round(algo.scores[[algo]]$error,2), "\nOffset: ", round(algo.scores[[algo]]$offset,2))) + 
      xlab("sd of cell types in bulks") + ylab("correlation") + geom_line() + ylim(0,1)
    plot(p1)
  }
  dev.off()
  
  # plot algorithm score vs correlation
  #cors <- apply(cor.mat, 2, mean)
  #cors <- cors[names(algo.scores)]
  
  #temp.df <- data.frame(algorithm = names(algo.scores), score = algo.scores, correlation = cors)
  #p2 <- ggplot(temp.df, aes(x = correlation, y = score)) + geom_line() + geom_point(aes(col = algorithm)) + ggtitle("Score vs. Correlation for each algorithm") + xlim(0,1)
  #pdf(paste0(plot.dir, "/score_vs_cor.pdf"), width = 16, height = 9)
  #plot(p2)
  #dev.off()
  
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
  if(!file.exists(paste0(temp.dir, "/results/deconv_results.RDS"))){
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
  }else{
    deconv.results <- readRDS(paste0(temp.dir, "/results/deconv_results.RDS"))
  }
  processed_results <- process_score_deconv_results(deconv.results, celltypes, temp.dir)
  
  results <- fit_model(processed_results, temp.dir, verbose)
  return(results$Scores)
}
