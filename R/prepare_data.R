prepare_data <- function(results.all, metric="cor") {
    df <- c()
    # extract real proportions from lists
    real.props <- results.all$bulk.props
    results.list <- results.all[-which(names(results.all) == "bulk.props")]
    if(length(results.list) == 1) results.list <- results.list[[1]]

    # iterate through the results list
    # depth of the list depends on the benchmark
    for(i in 1:length(results.list)) {
        results <- results.list[[i]]
        for(res in results) {
            # if the list has three levels the top level represents the geneset or amount of samples used
            # if res contains lists, it is either geneset or sample benchmark
            if(is.list(res[[1]])) {
                for(r in res) {
                    scores <- c()
                    name <- r$name
                    time <- as.numeric(r$times)

                    # try to determine what information is present in the lists
                    # this should do for now, but should be more general
                    # for further benchmarking modes
                    if(any(is.na(as.numeric(names(results.list)[i])))) {
                        geneset <- names(results.list)[i]
                        fraction <- 100
                    }else{
                        geneset <- NA
                        fraction <- as.numeric(names(results.list)[i])
                    }

                    # calculate condition number if reference is available
                    if(!is.null(r$sig.matrix)) {
                        cond.num <- kappa(r$sig.matrix, exact = F)
                    }else{
                        cond.num <- NA
                    }
		    
                    # if the deconvolution worked evaluate according to given metric
                    if(!all(is.null(r$est.props)) && !all(is.na(r$est.props))){
                        # performance per cell type
                        for (t in intersect(rownames(r$est.props), rownames(real.props))) {
                            temp.score <- evaluate_deconvolution(r$est.props[t,], real.props[t,])[[metric]]
                            scores <- c(scores, temp.score)
                            df <- rbind(df, c(name, temp.score, t, geneset, metric, time, fraction, cond.num))
                        }
                        # deal with NAs here...
                        if(any(is.na(scores))) {
                            if(metric == "cor")
                                scores[is.na(scores)] <- 0
                        }
                        # overall performance (average over per-cell-type-performance); store as cell type "overall"
                        df <- rbind(df, c(name, mean(scores), "overall", geneset, metric, time, fraction, cond.num))
                    }
                }
            }else{
                # this is executed if the list has two levels -> bulk benchmark
                # this is actually the same as above, only with fixed values for fraction and gene set
                # would be more elegant to join both cases but it works for now
                scores <- c()
                r <- res
                name <- r$name
                time <- as.numeric(r$times)
                if(!is.null(r$sig.matrix)) {
                        cond.num <- kappa(r$sig.matrix, exact = F)
                    }else{
                        cond.num <- NA
                    }
                if(!all(is.null(r$est.props)) && !all(is.na(r$est.props))){
                    for (t in intersect(rownames(r$est.props), rownames(real.props))) {
                        temp.score <- evaluate_deconvolution(r$est.props[t,], real.props[t,])[[metric]]
                        scores <- c(scores, temp.score)
                        df <- rbind(df, c(name, temp.score, t, NA, metric, time, 100, cond.num))
                    }
                    # deal with NAs here...
                    if(any(is.na(scores))) {
                            if(metric == "cor")
                                scores[is.na(scores)] <- 0
                        }
                    df <- rbind(df, c(name, mean(scores), "overall", NA, metric, time, 100, cond.num))
                }
            }
        }
    }
    df <- as.data.frame(df)
    if(ncol(df) != 8) 
	    return(data.frame())
    colnames(df) <- c("algorithm", "score", "cell_type", "geneset", "metric", "time", "fraction", "condition_number")
    df$score <- as.numeric(as.character(df$score))
    df$time <- as.numeric(as.character(df$time))
    df$condition_number <- as.numeric(as.character(df$condition_number))
    df$fraction <- as.numeric(as.character(df$fraction))
    return(df)
}
