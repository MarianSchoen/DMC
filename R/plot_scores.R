evaluation_plot <- function(results.df, title = NULL, metric = "cor", real.props = NULL, celltype.order = NULL, algorithm.order = NULL) {
    if(!is.data.frame(results.df)){
      stop("results.df must be a data frame")
    }
    if(!all(c("algorithm", "score", "cell_type", "metric") %in% colnames(results.df))){
      stop("required columns missing from results.df")
    }
    if (is.null(title)) {
        title <- "deconvolution quality"
    }
    if(!is.null(real.props) && !is.matrix(real.props)){
      stop("real.props is not a matrix")
    }
    # reduce to one entry per cell type and algorithm by taking the mean
    quality.scores <- c()
    for (a in unique(results.df$algorithm)) {
      for (ct in unique(results.df$cell_type)) {
        subset <- which(results.df$algorithm == a & results.df$cell_type == ct)
        if(length(subset > 0)) {
            quality.scores <- rbind(quality.scores, 
                                    c(mean(results.df[subset,]$score),
                                    sd(results.df[subset,]$score),
                                    ct, a, metric))
        }
      }
    }
    
    if (is.null(nrow(quality.scores))) {
        warning("No scores to plot")
        return(NULL)
    }
      # convert to data frame for later plotting
      quality.scores <- as.data.frame(quality.scores)
      colnames(quality.scores) <- c("value",
                                    "sd",
                                    "cell_type",
                                    "algorithm",
                                    "measure")
      
      # create labels for cell types in table plots
      if(is.null(celltype.order)){
        labels <- levels(quality.scores$cell_type)
        labels <- c("overall", labels[-which(labels == "overall")])
      }else{
        labels <- celltype.order
      }
      quality.scores$cell_type <- factor(quality.scores$cell_type, levels = labels)
      if(!is.null(real.props)) {
        for (i in 1:length(labels)) {
              if(labels[i] != "overall"){
                      prop <- round(mean(real.props[labels[i],]), 2)
                      labels[i] <-
                          paste(labels[i], "\n", prop * 100, "%", sep = "")
              }
            }
        }

      # cleaner plot by setting negative correlations to 0
      quality.scores$value <- as.numeric(as.character(quality.scores$value))
      quality.scores$value[quality.scores$value < 0] <- 0
      quality.scores$value[is.na(quality.scores$value)] <- 0

      # order by performance if ordering is not given
      ranking <- tapply(quality.scores$value, quality.scores$algorithm, mean)
      sds <- tapply(quality.scores$value, quality.scores$algorithm, sd)
      if(is.null(algorithm.order)){
        quality.scores$algorithm <- factor(quality.scores$algorithm,
                                           levels = levels(quality.scores$algorithm)[order(ranking)])
        sds <- sds[order(ranking)]
        ranking <- sort(ranking)
      }else{
        quality.scores$algorithm <- factor(quality.scores$algorithm, levels = algorithm.order)
        sds <- sds[algorithm.order]
        ranking <- ranking[algorithm.order]
      }
      

      # create table plot
      score.plot <- ggplot(quality.scores,
                         aes(x = as.numeric(cell_type), y = as.numeric(algorithm))) +
        geom_point(aes(
          size = value,
          col = value,
          fill = value
        ), shape = 22) +
        geom_point(
          shape = 22,
          fill = NA,
          color = "black",
          size = 25
        ) +
        xlab("cell type") + 
        ylab("algorithm") +
        ggtitle(title, subtitle = metric) +
        scale_size_continuous(
          limits = c(0, 1),
          range = c(min(quality.scores$value)*25, max(quality.scores$value) * 25),
          guide = "none"
        ) +
        theme(
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          title = element_text(size = 18),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)
        ) +
        scale_y_continuous(
          breaks = 1:length(levels(quality.scores$algorithm)),
          limits = c(0.5, length(levels(
            quality.scores$algorithm
          )) + 0.5),
          labels = paste(
            levels(quality.scores$algorithm),
            "\nr=",
            round(ranking, 2),
            " +/- ",
            round(sds, 2),
            sep = ""
          ),
          minor_breaks = seq(
            0,
            length(levels(
              quality.scores$algorithm
            )) + 1,
            0.1)
        ) +
        scale_x_continuous(
          breaks = 1:length(levels(quality.scores$cell_type)),
          limits = c(0.5, length(levels(
            quality.scores$cell_type
          )) + 0.5),
          labels = labels,
          minor_breaks = seq(
            0,
            length(levels(
              quality.scores$cell_type
            )) + 1,
            0.1)
        ) +
        scale_color_gradient2(
          low = "red",
          mid = "orange",
          high = "green",
          midpoint = 0.5
        ) +
        scale_fill_gradient2(
          low = "red",
          mid = "orange",
          high = "green",
          midpoint = 0.5
        )
    # return the plot and the ordering of the axes
    return(list(plot = score.plot, celltype.order = levels(quality.scores$cell_type), algorithm.order = levels(quality.scores$algorithm)))
}
