#' create condition number plots for deconvolution results
#' 
#' @param results.df data frame containing results as returned by prepare_data
#' @param metric evaluation metric; either string 'cor' (default) or a function
#' @param metric.name string, name of the evaluation metric used; not needed if metric is a string ("cor"). If metric is a function and metric.name
#' is NULL, the default will be "custom metric"
#' @param algorithm.order character vector of algorithm names specifying the plotting order
#' @return list containing 3 plots:
#' 1) cond_num_plot - barplot of average condition numbers
#' 2) cond_vs_score - scatter plot of average condition number vs. score for each algorithm
#' 3) variation_plot - standard deviation of score vs. standard deviation of condition number

plot_cond_num <- function(results.df, metric = "cor", metric.name = NULL, algorithm.order = NULL){
    if(!is.data.frame(results.df)){
        stop("results.df must be a data frame")
    }
    if(!all(c("score", "algorithm", "cell_type", "condition_number") %in% colnames(results.df))){
        stop("required columns missing from results.df")
    }
    if(is.character(metric)){
    if(metric != "cor"){
			stop("metric must be either \"cor\" or a function")
		}else{
            if(is.null(metric.name) || !is.character(metric.name)){
				metric.name <- "custom metric"
			}
			metric <- cor
		}
  }else{
    if(!is.function(metric)){
			stop("Function corresponding to 'metric' could not be found.")
		}else{
            if(is.null(metric.name) || !is.character(metric.name)){
				metric.name <- "custom metric"
			}
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
    # use only 'overall' rows
    results.df <- results.df[which(results.df$cell_type == "overall"), ]
    overall.df <- c()

    # create means and sds of scores and condition numbers from repetitions
    for(a in unique(results.df$algorithm)){
        overall.df <- rbind(overall.df, 
        c(
            a, 
            mean(results.df[which(results.df$algorithm == a),"condition_number"]), 
            sd(results.df[which(results.df$algorithm == a),"condition_number"]),
            mean(results.df[which(results.df$algorithm == a),"score"]),
            sd(results.df[which(results.df$algorithm == a),"score"])
            )
        )
    }
    overall.df <- as.data.frame(overall.df)
    colnames(overall.df) <- c("algorithm", 
    "condition_number", 
    "condition_variation", 
    "score" ,
    "score_variation"
    )

    # fix variable types
    overall.df$condition_number <- as.numeric(as.character(overall.df$condition_number))
    overall.df$condition_variation <- as.numeric(as.character(overall.df$condition_variation))
    overall.df$score <- as.numeric(as.character(overall.df$score))
    overall.df$score_variation <- as.numeric(as.character(overall.df$score_variation))

    # order the factor levels for plot order
    if(!is.null(algorithm.order)){
        overall.df$algorithm <- factor(overall.df$algorithm, levels = algorithm.order)
    }
    
    # create plot labels (for barplot)
    cond_labs <- levels(overall.df$algorithm)
    for(i in 1:length(cond_labs)){
        if(is.na(overall.df[which(overall.df$algorithm == cond_labs[i]), "condition_number"])){
            cond_labs[i] <- paste(cond_labs[i], "\n data missing", sep ="")
        }
    }

    # plot condition numbers as barplot
    cond_num_plot <- ggplot(overall.df) +
        geom_bar(aes(x = algorithm, y = condition_number, fill = algorithm), stat = "identity", position = "dodge") +
        geom_errorbar(aes(x = algorithm, ymin = condition_number - condition_variation, ymax = condition_number + condition_variation), width = 0.2) +
        ggtitle("average signature matrix condition number") +
        ylab("condition number") +
        xlab("algorithm") +
        scale_x_discrete(limits = levels(overall.df$algorithm), labels = cond_labs) +
        theme(
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 22),
            title = element_text(size = 24),
            axis.title.x = element_text(size = 22),
            axis.text.x = element_text(size = 20),
            axis.title.y = element_text(size = 22),
            axis.text.y = element_text(size = 20)
      )

    # plot score vs condition number
    cond_vs_score <- ggplot(overall.df) +
        geom_point(aes(x = condition_number, y = score, col = algorithm), size = 3, na.rm = T) +
        ggtitle("score vs condition number") +
        xlab("condition number") +
        ylab("score") + ylim(0,1) +
        theme(
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 22),
            title = element_text(size = 24),
            axis.title.x = element_text(size = 22),
            axis.text.x = element_text(size = 20),
            axis.title.y = element_text(size = 22),
            axis.text.y = element_text(size = 20)
      )

    # plot sd of score vs sd of condition number
    variation_plot <- ggplot(overall.df) +
        geom_point(aes(x = condition_variation, y = score_variation, col = algorithm), size = 3, na.rm = T) +
        ggtitle("SDs of score vs SDs of condition number") +
        xlab("condition number SD") +
        ylab("score SD") +
        theme(
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 22),
            title = element_text(size = 24),
            axis.title.x = element_text(size = 22),
            axis.text.x = element_text(size = 20),
            axis.title.y = element_text(size = 22),
            axis.text.y = element_text(size = 20)
      )

    return(list(cond_num_plot = cond_num_plot, cond_vs_score = cond_vs_score, variation_plot = variation_plot))
}
