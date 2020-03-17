plot_cond_num <- function(results.df){
    require(ggplot2)
    # do this only for 'overall' rows
    results.df <- results.df[which(results.df$cell_type == "overall"), ] 
    overall.df <- c()
    print(results.df)
    print(str(results.df))

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
    print(overall.df)
    print(str(overall.df))
    overall.df$condition_number <- as.numeric(as.character(overall.df$condition_number))
    overall.df$condition_variation <- as.numeric(as.character(overall.df$condition_variation))
    overall.df$score <- as.numeric(as.character(overall.df$score))
    overall.df$score_variation <- as.numeric(as.character(overall.df$score_variation))
    print(str(overall.df))
    cond_num_plot <- ggplot(overall.df) +
        geom_bar(aes(x = as.numeric(algorithm), y = condition_number, fill = algorithm), stat = "identity", position = "dodge") +
        ggtitle("average signature matrix condition number") +
        ylab("condition number") +
        xlab("algorithm")

    cond_vs_score <- ggplot(overall.df) +
        geom_point(aes(x = condition_number, y = score, col = algorithm)) +
        ggtitle("score vs condition number") +
        xlab("condition number") +
        ylab("score") +
        ylim(0,1)

    variation_plot <- ggplot(overall.df) +
        geom_point(aes(x = condition_variation, y = score_variation, col = algorithm)) +
        ggtitle("SDs of score vs SDs of condition number") +
        xlab("condition number SD") +
        ylab("score SD")

    return(list(cond_num_plot = cond_num_plot, cond_vs_score = cond_vs_score, variation_plot=variation_plot))
}