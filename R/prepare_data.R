prepare_data <- function(results.all, metric="cor") {
    df <- c()
    real.props <- results.all$bulk.props
    results.list <- results.all[-which(names(results.all) == "bulk.props")]
    if(length(results.list) == 1) results.list <- results.list[[1]]
    # iterate through the results list
    for(i in 1:length(results.list)) {
        results <- results.list[[i]]
        for(res in results) {
            # if the list has three levels the top level represents the geneset used
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
                        fraction = 100
                    }else{
                        geneset <- NA
                        fraction <- as.numeric(names(results.list)[i])
                    }

                    if(!is.null(r$sig.matrix)) {
                        cond.num <- kappa(r$sig.matrix, exact = F)
                    }else{
                        cond.num <- NA
                    }
		    
                    if(!is.null(r$est.props) && !is.na(r$est.props)){
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
                        df <- rbind(df, c(name, mean(scores), "overall", geneset, metric, time, fraction, cond.num))
                    }
                }
            }else{
                scores <- c()
                r <- res
                name <- r$name
                time <- as.numeric(r$times)
                if(!is.null(r$sig.matrix)) {
                        cond.num <- kappa(r$sig.matrix, exact = F)
                    }else{
                        cond.num <- NA
                    }
                if(!is.null(r$est.props) && !is.na(r$est.props)){
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
    colnames(df) <- c("algorithm", "score", "cell_type", "geneset", "metric", "time", "fraction", "condition_number")
    df$score <- as.numeric(as.character(df$score))
    df$time <- as.numeric(as.character(df$time))
    df$condition_number <- as.numeric(as.character(df$condition_number))
    return(df)
}
