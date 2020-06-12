#' create score plot
#' 
#' @param results.df data frame as returned by prepare_data
#' @param title character, plot title
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @param real.props non-negative numeric matrix, with cell types as rows, 
#' and bulk RNA-Seq profiles.
#' @param celltype.order character vector of cell types specifying the plotting order
#' @param algorithm.order character vector of algorithm names specifying the plotting order (left to right)
#' @result list with 3 entries:
#' 1) plot - table plot of deconvolution results
#' 2) celltype.order - character vector containing cell types in plotting order
#' 3) algorithm.order - character vector containing algorithms in plotting order

evaluation_plot <- function(results.df, title = NULL, metric = "cor", real.props = NULL, celltype.order = NULL, algorithm.order = NULL) {
  # parameter check
    if(!is.data.frame(results.df)){
      stop("results.df must be a data frame")
    }
    if(!all(c("algorithm", "score", "cell_type", "metric") %in% colnames(results.df))){
      stop("required columns missing from results.df")
    }
    if (is.null(title)) {
        title <- "deconvolution quality"
    }
    if(is.character(metric)){
    if(metric != "cor"){
			stop("metric must be either \"cor\" or a function")
		}else{
			metric <- cor
		}
  }else{
    if(!is.function(metric)){
			stop("Function corresponding to 'metric' could not be found.")
		}
  }
    if(!is.null(real.props) && !is.matrix(real.props)){
      stop("real.props is not a matrix")
    }
    if(!is.null(celltype.order)){
        if(!is.character(celltype.order)){
            stop("celltype.order must be a charcter vector")
        }
        if(!all(celltype.order %in% unique(results.df$cell_type)) || length(celltype.order) != length(unique(results.df$cell_type))){
            stop("celltype.order does not fit the cell_type column of results.df")
        }
    }
    if(!is.null(algorithm.order)){
        if(!is.character(algorithm.order)){
            stop("celltype.order must be a charcter vector")
        }
        if(!all(algorithm.order %in% unique(results.df$algorithm)) || length(algorithm.order) != length(unique(results.df$algorithm))){
            stop("algorithm.order does not fit the algorithm column of results.df")
        }
    }

    # reduce to one entry per cell type and algorithm by taking the mean over repetitions
    quality.scores <- c()
    for (a in unique(results.df$algorithm)) {
      for (ct in unique(results.df$cell_type)) {
        subset <- which(results.df$algorithm == a & results.df$cell_type == ct)
        if(length(subset > 0)) {
            quality.scores <- rbind(quality.scores, 
                                    c(mean(results.df[subset,]$score),
                                    sd(results.df[subset,]$score),
                                    ct, a, as.character(substitute(metric))))
        }
      }
    }
    
    if (is.null(nrow(quality.scores))) {
        warning("No scores to plot")
        return(NULL)
    }
    # convert to data frame for plotting
    quality.scores <- as.data.frame(quality.scores)
    colnames(quality.scores) <- c("value",
                                  "sd",
                                  "cell_type",
                                  "algorithm",
                                  "measure")
    # create labels for cell types in table plots
    if(is.null(celltype.order)){
      labels <- unique(as.character(quality.scores$cell_type))
      # make sure 'overall' is the first type
      labels <- c("overall", labels[-which(labels == "overall")])
    }else{
      labels <- celltype.order
    }
    quality.scores$cell_type <- factor(quality.scores$cell_type, levels = labels)

    # if real.props is available, add average fraction of this cell type in the bulks to cell type label
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
    quality.scores$sd <- as.numeric(as.character(quality.scores$sd))
    if(any(quality.scores$value < 0)){
      quality.scores$value[quality.scores$value < 0] <- 0
    }
    if(any(is.na(quality.scores$value))){
      warning("NAs detected")
      quality.scores$value[is.na(quality.scores$value)] <- 0
    }

    # order algorithms by performance if ordering is not given
    algos <- quality.scores$algorithm[which(quality.scores$cell_type == "overall")]
    ranking <- quality.scores$value[which(quality.scores$cell_type == "overall")]
    sds <- quality.scores$sd[which(quality.scores$cell_type == "overall")]
    names(ranking) <- algos
    names(sds) <- algos

    if(is.null(algorithm.order)){
      quality.scores$algorithm <- factor(quality.scores$algorithm,
                                          levels = names(ranking)[order(ranking)])
      sds <- sds[order(ranking)]
      ranking <- sort(ranking)
    }else{
      # order by algorithm.order
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
      ggtitle(title, subtitle = as.character(substitute(metric))) +
      scale_size_continuous(
        limits = c(0, 1),
        #range = c(min(quality.scores$value)*25, max(quality.scores$value) * 25),
        range = c(0,25),
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
          0.5,
          length(levels(
            quality.scores$algorithm
          )) + 0.5,
          0.1)
      ) +
      scale_x_continuous(
        breaks = 1:length(levels(quality.scores$cell_type)),
        limits = c(0.5, length(levels(
          quality.scores$cell_type
        )) + 0.5),
        labels = labels,
        minor_breaks = seq(
          0.5,
          length(levels(
            quality.scores$cell_type
          )) + 0.5,
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
