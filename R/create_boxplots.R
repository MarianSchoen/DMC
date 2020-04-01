# create boxplots from a data frame as returned by prepare_data()
create_boxplots <- function(results.df, metric = "cor", celltype.order = NULL, algorithm.order = NULL) {
    require(ggplot2)
    if(!is.data.frame(results.df)){
        stop("results.df must be a data frame")
    }
    if(!all(c("score", "algorithm", "cell_type") %in% colnames(results.df))){
        stop("required columns missing from results.df")
    }
    if(!metric %in% c("cor", "mad", "rmsd")){
        stop("unknown metric. choose one of 'cor', 'mad', 'rmsd'")
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
    overall.df <- results.df[which(results.df$cell_type == "overall"), ]
    
    # order algorithms by performance or by given vector
    if(is.null(algorithm.order)){
        performances <- tapply(overall.df$score, overall.df$algorithm, median)
        results.df$algorithm <- factor(results.df$algorithm, levels = levels(overall.df$algorithm)[order(performances)])
    }else{
        results.df$algorithm <- factor(results.df$algorithm, levels = algorithm.order) 
    }
    # order cell types by supplied vector
    if(!is.null(celltype.order)){
        results.df$cell_type <- factor(results.df$cell_type, levels = celltype.order)
    }
    
    
    metric <- results.df$metric[1]

    # create boxplots for each cell type (including overall)
    cell.type.plots <- list()
    for(t in levels(results.df$cell_type)) {
        sub.df <- results.df[which(results.df$cell_type == t), ]
        cell.type.plots[[t]] <- ggplot(sub.df, aes(x = algorithm, y = score)) +
            geom_boxplot(aes(col = algorithm)) +
            xlab("algorithm") +
            ylab(metric) +
            ggtitle("quality of deconvolution results", subtitle = t)+
            theme(
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 16),
                title = element_text(size = 22),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 20),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)
            ) + 
            scale_x_discrete(limits = levels(results.df$algorithm))
        if(metric == "cor"){
            cell.type.plots[[t]] <- cell.type.plots[[t]] + ylim(0,1)
        }
        if(t == "overall") {
            cell.type.plots[[t]] <- cell.type.plots[[t]] + ylab(paste("average overall", metric))
        }
    }
    return(list(cell.type.plots = cell.type.plots))
}