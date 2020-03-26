plot_cond_num <- function(results.df, metric = "cor"){
    require(ggplot2)
    if(!is.data.frame(results.df)){
        stop("results.df must be a data frame")
    }
    if(!all(c("score", "algorithm", "cell_type", "condition_number") %in% colnames(results.df))){
        stop("required columns missing from results.df")
    }
    if(!metric %in% c("cor", "mad", "rmsd")){
        stop("unknown metric. must be one of 'cor', 'mad', 'rmsd'")
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

    # plot condition numbers as barplot
    cond_num_plot <- ggplot(overall.df) +
        geom_bar(aes(x = as.numeric(algorithm), y = condition_number, fill = algorithm), stat = "identity", position = "dodge") +
        ggtitle("average signature matrix condition number") +
        ylab("condition number") +
        xlab("algorithm")

    # plot score vs condition number
    cond_vs_score <- ggplot(overall.df) +
        geom_point(aes(x = condition_number, y = score, col = algorithm)) +
        ggtitle("score vs condition number") +
        xlab("condition number") +
        ylab("score")
    if(metric == "cor"){
        cond_vs_score <- cond_vs_score + ylim(0,1)
    }

    # plot sd of score vs sd of condition number
    variation_plot <- ggplot(overall.df) +
        geom_point(aes(x = condition_variation, y = score_variation, col = algorithm)) +
        ggtitle("SDs of score vs SDs of condition number") +
        xlab("condition number SD") +
        ylab("score SD")

    return(list(cond_num_plot = cond_num_plot, cond_vs_score = cond_vs_score, variation_plot=variation_plot))
}