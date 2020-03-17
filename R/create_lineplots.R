create_lineplots <- function(results.df, genesets = NULL, available.features = NULL) {
    require(ggplot2)
    overall.df <- results.df[which(results.df$cell_type == "overall"), ]
    # order algorithms by performance
    performances <- tapply(overall.df$score, overall.df$algorithm, median)
    overall.df$algorithm <- factor(overall.df$algorithm, levels = levels(overall.df$algorithm)[order(performances)])
    results.df$algorithm <- factor(results.df$algorithm, levels = levels(overall.df$algorithm)[order(performances)])

    metric <- overall.df$metric[1]
    if(!is.null(genesets)) {
        if(all(results.df$geneset %in% names(genesets))){
            geneset.labs <- paste(names(genesets), "\n(", as.numeric(sapply(genesets, function(x) length(which(x %in% available.features)))), " genes)", sep = "")
        }else{
            geneset.labs <- levels(results.df$geneset)
        }
        geneset.limits <- names(sort(sapply(genesets, function(x) length(which(x %in% available.features)))))
    }else{
        geneset.labs <- levels(results.df$geneset)
        geneset.limits <- levels(results.df$geneset)
    }
   
    
    # create overall performance plot
    #correlation.plot <- ggplot(overall.df, aes(
    #x = geneset,
    #y = score,
    #group = algorithm, col = algorithm
    #)) +
    #geom_line(size = 2) + geom_point() +
    #xlab("gene set") + ylab(paste("average",metric)) +
    #ggtitle("deconvolution quality using different gene sets", subtitle = metric) +
    #ylim(0, 1) +
    #theme(
    #    legend.text = element_text(size = 20),
    #    legend.title = element_text(size = 22),
    #    title = element_text(size = 24),
    #    axis.title.x = element_text(size = 24),
    #    axis.text.x = element_text(size = 20),
    #    axis.title.y = element_text(size = 24),
    #    axis.text.y = element_text(size = 20)
    #) +
    #scale_x_discrete(limits = geneset.limits, labels = geneset.labs) +
    #guides(linetype = guide_legend(override.aes = list(size = 2)))

    # plot runtime
    
    # create plot per cell type
    cell.type.plots <- list()
    for(t in unique(results.df$cell_type)) {
        #if(t == "overall") next
        sub.df <- results.df[which(results.df$cell_type == t), ]
	temp.scores <- tapply(sub.df$score, list(sub.df$algorithm, sub.df$geneset), mean)
    	temp.times <- tapply(sub.df$time, list(sub.df$algorithm, sub.df$geneset), mean)
    	sub.df <- c()
    	for(i in 1:ncol(temp.scores)) {
	    for(j in 1:nrow(temp.scores)) {
		    sub.df <- rbind(sub.df,
		    			c(rownames(temp.scores)[j], colnames(temp.scores)[i], temp.scores[j, i], temp.times[j,i])
		    )
	    }
    	}
	
    	sub.df <- as.data.frame(sub.df)
    	colnames(sub.df) <- c("algorithm", "geneset", "score", "time")
    	sub.df$score <- as.numeric(as.character(sub.df$score))
    	sub.df$time <- as.numeric(as.character(sub.df$time))
        cell.type.plots[[t]] <- ggplot(sub.df, aes(x=geneset, y = score, group = algorithm, col = algorithm)) +
            geom_line(size = 2) + geom_point() +
            xlab("gene set (increasing size)") +
            ylab(metric) +
            ggtitle(paste(
            "deconvolution quality using different gene sets (", t, ")",
            sep = ""
            )) +
            ylim(0, 1) +
            scale_x_discrete(limits = geneset.limits) +
            guides(linetype = guide_legend(override.aes = list(size = 2)))
    	if(t == "overall") {
		runtime.plot <- ggplot(sub.df, aes(
      		x = geneset, y = log(time, 10),
      		group = algorithm, col = algorithm
    		)) +
      		geom_line(size = 2) + geom_point() +
      		xlab("gene set") + ylab("log time (s)") +
      		ggtitle("runtime of algorithms using different gene sets") +
      		theme(
        		legend.text = element_text(size = 20),
        		legend.title = element_text(size = 22),
        		title = element_text(size = 24),
        		axis.title.x = element_text(size = 22),
        		axis.text.x = element_text(size = 20),
        		axis.title.y = element_text(size = 22),
        		axis.text.y = element_text(size = 20)
      		) +
      		scale_x_discrete(limits = geneset.limits, labels = geneset.labs) +
      		guides(linetype = guide_legend(override.aes = list(size = 2)))
	}
    }
    return(list(runtime.plot = runtime.plot, cell.type.plots = cell.type.plots))
}