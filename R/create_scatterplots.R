# written by Tim Mirus

create_scatterplots <- function(results.list, real.props = NULL, training.pheno = NULL, real = FALSE) {
  require(ggplot2)
  real.props <- results.list$bulk.props
    # create scatter plots of real vs estimate props
    if (!is.null(real.props) && !is.null(training.pheno) && real) {
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
    # create scatter plots only for one repetition
    for (res in results.list[[1]]) {   
      # only cell types that are in real props and estimates
      cts <- intersect(rownames(res$est.props), rownames(real.props))
      # create data frame for plotting (containing real and estimates for each cell type)
      df <- c()
      if(!length(cts) > 0) next
      for (i in 1:length(cts)) {
        t <- cts[i]
        temp.df <- data.frame(real = real.props[t, ], estimate = res$est.props[t, ])
        rownames(temp.df) <- colnames(real.props)
        df <- rbind(df, cbind(temp.df, t))
      }
      df <- data.frame(df)
      colnames(df) <- c('real', 'estimate', 'type')
      df$real <- as.numeric(as.character(df$real))
      df$estimate <- as.numeric(as.character(df$estimate))
      
      # add pearson r to cell type titles
      labs <- levels(df$type)
      temp <- rep("", length(labs))
      for (i in 1:length(labs)) {
        data.points <- which(df$type == labs[i])
        temp[i] <- paste(labs[i], "\nr = ", round(cor(df[data.points, 1], df[data.points, 2]), 2), sep = "")
      }
      names(temp) <- labs
      labs <- temp
      
      # create scatter plots
      scatter.plot <- ggplot(df) +
        geom_point(aes(x = real, y = estimate, col = type), size = 2) +
        facet_grid(
          rows = NULL,
          cols = vars(type),
          labeller = labeller(type = labs)
        ) +
        xlab("real") + ylab("estimate") +
        labs(col = 'cell type') +
        theme(
          strip.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)
        ) + ggtitle(res$name)
      scatter.plots[[res$name]] <- scatter.plot
    }
      return(scatter.plots)
}
