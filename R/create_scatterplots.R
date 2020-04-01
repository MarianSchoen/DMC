create_scatterplots <- function(results.list, real.props = NULL, training.pheno = NULL, real = FALSE, celltype.order = NULL, algorithm.order = NULL) {
  require(ggplot2)
  if(!is.list(results.list)){
    stop("results.list must be a list as returned by deconvolute()")
  }
  if(is.null(real.props)){
    if(!"bulk.props" %in% names(results.list)){
      stop("results.list must contain an entry bulk.props; alternatively supply real.props")
    }
    real.props <- results.list$bulk.props
  }
  
    # create scatter plots of real vs estimate props
    if (!is.null(real.props) && !is.null(training.pheno) && real && "patient" %in% colnames(training.pheno)) {
      bulk.shapes <- as.factor(ifelse(
        colnames(real.props) %in% training.pheno$patient,
        "training",
        "test only"
      ))
    } else {
      bulk.shapes <- rep(1, ncol(real.props))
    }

    results.list <- results.list[["results.list"]]
    scatter.plots <- list()

    dfs <- list()
    for(n in names(results.list[[1]])){
      dfs[[n]] <- data.frame()
    }
    # create scatter plots only for one repetition
    for (j in 1:length(results.list)) {   
      for(res in results.list[[j]]){
        # only cell types that are in real props and estimates
        cts <- intersect(rownames(res$est.props), rownames(real.props))
        # create data frame for plotting (containing real and estimates for each cell type)
        df <- c()
        if(!length(cts) > 0) next
        for (i in 1:length(cts)) {
          t <- cts[i]
          temp.df <- data.frame(real = real.props[t, ], estimate = res$est.props[t, ])
          rownames(temp.df) <- colnames(real.props)
          df <- rbind(df, cbind(temp.df, t, j, evaluate_deconvolution(real.props[t, ], res$est.props[t, ])$cor))
        }
        df <- data.frame(df)
        colnames(df) <- c('real', 'estimate', 'type', 'repetition', 'cor')
        df$real <- as.numeric(as.character(df$real))
        df$estimate <- as.numeric(as.character(df$estimate))
        dfs[[res$name]] <- rbind(dfs[[res$name]], df)
      }
    }
    # plot here...
    for(name in names(dfs)){
      df <- dfs[[name]]
      if(!is.null(celltype.order)){
          df$type <- factor(df$type, levels = celltype.order)
        }
      cors <- as.numeric(as.character(unique(df$cor)))
      
      labs <- levels(df$type)
      temp <- rep("", length(labs))
      for(i in 1:length(labs)){
        data.points <- which(df$type == labs[i])
        temp[i] <- paste(labs[i], "\nr = ", round(mean(as.numeric(as.character(unique(df[data.points, "cor"])))), 2), " +/- ", round(sd(as.numeric(as.character(unique(df[data.points, "cor"])))), 2), sep = "")
      }
      names(temp) <- labs
      labs <- temp

      scatter.plot <- ggplot(df) +
        geom_point(aes(x = real, y = estimate, col = as.factor(repetition)), size = 2) +
        facet_grid(
          rows = NULL,
          cols = vars(type),
          labeller = labeller(type = labs)
        ) +
        xlab("real") + ylab("estimate") +
        labs(col = 'repetition') +
        theme(
            strip.text.x = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16)
        ) + 
        ggtitle(res$name)
        scatter.plots[[name]] <- scatter.plot
    }

     # # add pearson r to cell type titles
        # labs <- levels(df$type)
        # temp <- rep("", length(labs))
        # for (i in 1:length(labs)) {
        #   data.points <- which(df$type == labs[i])
        #   temp[i] <- paste(labs[i], "\nr = ", round(cor(df[data.points, 1], df[data.points, 2]), 2), sep = "")
        # }
        # names(temp) <- labs
        # labs <- temp
        
        # # create scatter plots
        # scatter.plot <- ggplot(df) +
        #   geom_point(aes(x = real, y = estimate, col = type), size = 2) +
        #   facet_grid(
        #     rows = NULL,
        #     cols = vars(type),
        #     labeller = labeller(type = labs)
        #   ) +
        #   xlab("real") + ylab("estimate") +
        #   labs(col = 'cell type') +
        #   theme(
        #     strip.text.x = element_text(size = 16),
        #     axis.title.x = element_text(size = 18),
        #     axis.title.y = element_text(size = 18),
        #     legend.title = element_text(size = 18),
        #     legend.text = element_text(size = 16)
        #   ) + ggtitle(res$name)
        # scatter.plots[[res$name]] <- scatter.plot
    if(!is.null(algorithm.order)){
      if(all(algorithm.order %in% names(scatter.plots)) && length(scatter.plots) == length(algorithm.order)){
        scatter.plots <- scatter.plots[algorithm.order]
      }
    }
      return(scatter.plots)
}
