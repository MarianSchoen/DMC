#' create runtime plot
#' 
#' @param results.df data frame containing deconvolution results as returned by prepare_data
#' @param title character, plot title
#' @param algorithm.order character vector of algorithm names specifying the plotting order
#' @return runtime plot

plot_runtime <- function(results.df, title = NULL, algorithm.order = NULL) {
    # parameter checks
    if(!is.data.frame(results.df)) {
        stop("Invalid input: results.df has to be data frame as returned by prepare_data()")
    }
    if(!all(c("algorithm", "time", "cell_type") %in% colnames(results.df))){
        stop("required columns missing in results.df")
    }
    if (is.null(title))
       title <- "runtime comparison"
    if(!is.null(algorithm.order)){
        if(!is.character(algorithm.order)){
            stop("celltype.order must be a charcter vector")
        }
        if(!all(algorithm.order %in% unique(results.df$algorithm)) || length(algorithm.order) != length(unique(results.df$algorithm))){
            stop("algorithm.order does not fit the algorithm column of results.df")
        }
    }
    # reduce to the 'overall' rows
    results.df <- results.df[which(results.df$cell_type == "overall"),]

    # average over all repetitions
    mean.times <- tapply(results.df$time, results.df$algorithm, mean)
    sd.times <- tapply(results.df$time, results.df$algorithm, sd)

    # create data frame containing runtimes
    runtimes <- data.frame(algorithm = names(mean.times), runtime = mean.times, sd = sd.times)
    if(!is.null(algorithm.order)){
        runtimes$algorithm <- factor(runtimes$algorithm, levels = algorithm.order)
    }
    
    # plot runtime as lineplot
    runtime.plot <- ggplot(runtimes) +
    geom_bar(
        aes(
        y = algorithm,
        x = log10(as.numeric(as.character(runtime)) + 1),
        fill = algorithm
        ),
        stat = "identity",
        position = "dodge"
    ) +
    ylab("algorithm") + xlab("time (s)") + ggtitle(title) +
    theme(
        legend.position = "none",
        title = element_text(size = 24),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20)
      ) +
        scale_x_continuous(
            labels = function(x){round(10^as.numeric(x) - 1)},
        )
    return(runtime.plot)
}