plot_runtime <- function(results.df, title = NULL) {
    require(ggplot2)
    if(!is.data.frame(results.df)) {
        stop("Invalid input: results.df has to be data frame as returned by prepare_data()")
    }
    if (is.null(title))
       title <- "runtime comparison"
    # reduce to the overall rows
    results.df <- results.df[which(results.df$cell_type == "overall"),]
    # average over all runs
    mean.times <- tapply(results.df$time, results.df$algorithm, mean)
    sd.times <- tapply(results.df$time, results.df$algorithm, sd)

    # create data frame and plot
    runtimes <- data.frame(algorithm = names(mean.times), runtime = mean.times, sd = sd.times)
    runtime.plot <- ggplot(runtimes) +
    geom_bar(
        aes(
        x = algorithm,
        y = log10(as.numeric(as.character(runtime)) + 1),
        fill = algorithm
        ),
        stat = "identity",
        position = "dodge"
    ) +
    xlab("algorithm") + ylab("time (logarithmic)") + ggtitle(title) +
    theme(
        legend.position = "None",
        title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)
    )
    return(runtime.plot)
}