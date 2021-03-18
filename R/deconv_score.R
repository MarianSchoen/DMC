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
      
      for(ct in celltypes){
	      if(ct %in% rownames(results[[i]]$real.props)){
        	ct.sd <- sd(results[[i]]$real.props[ct,])
	      }else{
		      ct.sd <- 0
	      }
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
      selection <- intersect(
        which(correlations$celltype == ct),
        which(correlations$algorithm == algo)
      )
      # make variance scale of cell types comparable
      if(! max(correlations[selection, "ct_sd"]) == 0){
        correlations[selection, "ct_sd"] <- (correlations[selection, "ct_sd"] / max(correlations[selection, "ct_sd"])) * 100
      }
      for(i in unique(correlations$index)){
        subselection <- intersect(
          selection,
          which(correlations$index == i)
        )
        m <- mean(correlations[subselection,"score"])
        e <- sd(correlations[subselection, "score"])
        ct_sd <- mean(correlations[subselection, "ct_sd"])
        #entropy <- mean(correlations[subselection, "entropy"])
        
        reduced.correlations <- rbind(
          reduced.correlations,
          data.frame(algorithm = algo, index = i, celltype = ct, score = m, sd = e, ct_sd = ct_sd)#, entropy = entropy)
        )
      }
    }
  }
  correlations <- reduced.correlations
  rm(reduced.correlations)
  saveRDS(correlations, paste0(results.dir, "/correlations_sd.RDS"))

  overall.fit.df <- data.frame()
  pdf(paste0(plot.dir, "/model_fits.pdf"), width = 8 * length(unique(correlations$algorithm)), height = 8)
  fits <- list()
  
  # 
  

    #fit.df <- data.frame() # contains results
    algo.res <- c() # algorithm score as string
    algo.names <- c() # algorithm names for algo.res
    for(algo in unique(correlations$algorithm)){
      if(verbose) cat(algo, "...\t")
      temp.cors <- correlations[which(correlations$algorithm == algo),]
      
      # get parameters of the model
      #
      # initial values
      params <- c(x_a = 1)
      for(ct in unique(correlations$celltype)){
        param_names <- c(names(params), paste0(ct, "_x"), paste0(ct, "_a"))
        params <- c(params, 0.1, 0.1)
        names(params) <- param_names
      }
      # function to optimize
      cost_fun <- function(p, df = temp.cors){
        residual <- 0
        for(ct in unique(df$celltype)){
          sub_df <- df[df$celltype == ct,]
          
          # weights for fit
          if(any(is.na(sub_df$sd))){
            w <- rep(1, nrow(sub_df))
          }else{
            w <- sub_df$sd
          }
          if(any(w == 0)){
            w[w==0] <- 1
          }
          
          # add residual for celltype ct
          residual <- residual + 
            sum(
                ((
                  sub_df$score - 
                    (
                      (1 / sqrt(1 + ((abs(p["x_a"]) + abs(p[paste0(ct, "_x")]))/sub_df$ct_sd)^2)) - 
                        abs(p[paste0(ct, "_a")])
                      )
                )/w)^2
              )
        }
        return(residual)
      }
      # optimize the cost function
      param_out <- optim(params, cost_fun, method = "BFGS")
      param_out$par <- abs(param_out$par)
      # print(param_out$par)
      
      fits[[algo]] <- param_out
      
      overall_x <- param_out$par["x_a"]
      fit.df <- data.frame()
      for(ct in unique(temp.cors$celltype)){
        ct_x <- param_out$par[paste0(ct, "_x")]
        ct_a <- param_out$par[paste0(ct, "_a")]
        sub_df <- temp.cors[which(temp.cors$celltype == ct),]
        
        predictions <- 1/sqrt(1 + ((ct_x+overall_x)/sub_df$ct_sd)^2) - ct_a
        
        fit.df <- rbind(fit.df,
                        data.frame(
                          correlation = sub_df$score,
                          sd = sub_df$ct_sd,
                          celltype = ct,
                          prediction = predictions,
                          algorithm = algo,
                          ct_error = ct_x,
                          offset = ct_a,
                          algo_error = overall_x,
                          score_sd = sub_df$sd
                        )
                        )
      }
      
      algo.names <- c(algo.names, algo)

      # error and offset
      algo.res <- c(algo.res, paste0("Algorithm SD(error): ", round(mean(fit.df$algo_error),4), "\nmean Cell Type SD(error): ", round(mean(fit.df$ct_error),4), "\nmean Offset: ", round(mean(fit.df$offset), 3)))
      overall.fit.df <- rbind(overall.fit.df, fit.df)
    }
    names(algo.res) <- algo.names
    overall.fit.df$algorithm <- factor(overall.fit.df$algorithm, levels = names(algo.res))
    
    if(nrow(overall.fit.df) > 0){
      # plot
      labs = paste(names(algo.res), algo.res, sep = "\n")
      names(labs) <- names(algo.res)
      for(algo in unique(overall.fit.df$algorithm)){
        sub_df <- overall.fit.df[which(overall.fit.df$algorithm == algo),]
        
        # removed a x=sd in geom_line here...
        p <- ggplot(sub_df, aes(x = sd, y = correlation)) + 
          geom_pointrange(aes(ymin = correlation - score_sd, ymax = correlation + score_sd), col = "black") + 
          geom_line(aes(y = prediction), col = "red") + 
          ggtitle("calculated correlation and fitted model", subtitle = labs[algo]) + 
          xlab("bulk data SD") + ylab("Correlation") + ylim(0,1)+
          facet_grid(cols = vars(celltype))
        plot(p)
      }
    }
    if(verbose) cat("\n")
###############################################
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
                          if(any(is.na(temp.df$ct_error))){
                            mean.sd <- mean(unique(temp.df$ct_error[!is.na(temp.df$ct_error)]))
                          }else{
                            mean.sd <- mean(unique(temp.df$ct_error))
                          }
                          
                          if(any(is.na(temp.df$offset))){
                            mean.offset <- mean(unique(temp.df$offset[!is.na(temp.df$offset)]))
                          }else{
                            mean.offset <- mean(unique(temp.df$offset))
                          }
                          
                          if(any(is.na(temp.df$algo_error))){
                            mean.error <- mean(unique(temp.df$algo_error[!is.na(temp.df$algo_error)]))
                          }else{
                            mean.error <- mean(unique(temp.df$algo_error))
                          }
                          return(list(error = mean.sd, offset = mean.offset, algo_error = mean.error))
                        }
  )
  names(algo.scores) <- unique(correlations$algorithm)

  if(verbose) {
    cat("Scores:\n")
    cat("\t\tCT_Error\tAlgo Error\tOffset\n")
    for(algo in names(algo.scores)){
      cat(algo, ": \t", algo.scores[[algo]]$error,"\t", algo.scores[[algo]]$algo_error,"\t", algo.scores[[algo]]$offset, "\n")
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
      ggtitle(paste0("Overall Correlation for ", algo), subtitle = paste0("Score [mean SD(Error)]: ", round(algo.scores[[algo]]$error,2), "\nScore 2 [mean Algorithm SD(Error)]: ", round(algo.scores[[algo]]$algo_error,2), "\nOffset: ", round(algo.scores[[algo]]$offset,2))) + 
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
  processed_results <- process_score_deconv_results(deconv.results, predict.types, temp.dir)
  
  results <- fit_model(processed_results, temp.dir, verbose)
  return(results$Scores)
}
