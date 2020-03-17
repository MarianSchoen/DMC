create_boxplots <- function(results.df, metric = "cor") {
    require(ggplot2)
    overall.df <- results.df[which(results.df$cell_type == "overall"), ]
    # order algorithms by performance
    performances <- tapply(overall.df$score, overall.df$algorithm, median)
    overall.df$algorithm <- factor(overall.df$algorithm, levels = levels(overall.df$algorithm)[order(performances)])
    results.df$algorithm <- factor(results.df$algorithm, levels = levels(overall.df$algorithm)[order(performances)])
    
    metric <- overall.df$metric[1]
    # plot results in boxplot
    overall.plot <- ggplot(overall.df, aes(x = algorithm, y = score)) +
    geom_boxplot(aes(col = algorithm)) +
    ylim(0, 1) +
    xlab("algorithm") +
    ylab(paste("average overall", metric)) +
    ggtitle("quality of deconvolution results",
        subtitle = "varying number of training samples"
    ) +
    theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)
    ) + 
    scale_x_discrete(limits = levels(overall.df$algorithm))

    cell.type.plots <- list()
    for(t in unique(results.df$cell_type)) {
        if(t == "overall") next
        sub.df <- results.df[which(results.df$cell_type == t), ]
	sub.df$algorithm <- factor(sub.df$algorithm, levels = levels(overall.df$algorithm))
        cell.type.plots[[t]] <- ggplot(sub.df, aes(x = algorithm, y = score)) +
            geom_boxplot(aes(col = algorithm)) +
            xlab("algorithm") +
            ylab(metric) +
            ggtitle("quality of deconvolution results", subtitle = t)
        if(metric == "cor"){
            cell.type.plots[[t]] <- cell.type.plots[[t]] + ylim(0,1)
        }
    }
    return(list(score.plot = overall.plot, cell.type.plots = cell.type.plots))
}
