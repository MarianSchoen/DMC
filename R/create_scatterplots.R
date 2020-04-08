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
    for(name in names(dfs)){
      df <- dfs[[name]]
      if(nrow(df) == 0) {
        scatter.plots[[name]] <- NULL
        next
      }
      if(!is.null(celltype.order)){
          df$type <- factor(df$type, levels = celltype.order)
        }
      cors <- matrix(NA, nrow = length(unique(df$repetition)), ncol = length(levels(df$type)))
      rownames(cors) <- 1:nrow(cors)
      colnames(cors) <- levels(df$type)
      for(r in unique(df$repetition)){
        for(ct in unique(df$type)){
          temp <- as.numeric(as.character(unique(df[which(df$repetition == r & df$type == ct),"cor"])))
          if(!is.null(temp) && !is.na(temp)){
            cors[r, ct] <- temp
          }
        }
      }
      
      labs <- levels(df$type)
      temp <- rep("", length(labs))
      for(i in 1:length(labs)){
        temp[i] <- paste(labs[i], "\nr = ", round(mean(cors[,labs[i]]), 2), " +/- ", round(sd(cors[,labs[i]]), 2), sep = "")
      }
      names(temp) <- labs
      labs <- temp
      
      scatter.plot <- ggplot(df, aes(x = real, y = estimate, col = as.factor(type), group = real)) +
        geom_boxplot(position = "dodge", varwidth = T) +
        geom_point(alpha = 0.1) +
        facet_grid(
          rows = NULL,
          cols = vars(type),
          labeller = labeller(type = labs)
        ) +
        xlab("real") + ylab("estimate") +
        labs(col = 'cell type') +
        theme(
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22),
          title = element_text(size = 24),
          axis.title.x = element_text(size = 22),
          axis.text.x = element_text(size = 20),
          axis.title.y = element_text(size = 22),
          axis.text.y = element_text(size = 20)
        ) + 
        ggtitle(name)
        scatter.plots[[name]] <- scatter.plot
    }
    if(!is.null(algorithm.order)){
      if(all(algorithm.order %in% names(scatter.plots)) && length(scatter.plots) == length(algorithm.order)){
        scatter.plots <- scatter.plots[algorithm.order]
      }
    }
      return(scatter.plots)
}
