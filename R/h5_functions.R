# written by Tim Mirus
write_data <- function(sc.counts = NULL, sc.pheno = NULL, bulk.counts = NULL, bulk.props = NULL, filename) {
    require(rhdf5)
    h5createFile(filename)
    # write singlecell stuff
    # assume that sc.counts and sc.pheno are never written independently...
    # in the context of this benchmark this is sensible
    if(!is.null(sc.counts) && !is.null(sc.pheno)){
        h5createGroup(filename, "singlecell")
        h5write(as.vector(rownames(sc.counts)), filename, "singlecell/geneids", write.attributes = TRUE)
        h5write(as.vector(colnames(sc.counts)), filename, "singlecell/cellids", write.attributes = TRUE)
        h5write(as.matrix(sc.counts), filename, "singlecell/data", write.attributes = TRUE)
        h5createGroup(filename, "singlecell/pheno")
        for(name in colnames(sc.pheno)){
            h5write(as.vector(sc.pheno[[name]]), filename, paste("singlecell/pheno", name, sep = "/"))
        }
    }
    # write bulks if present
    if(!is.null(bulk.counts)){
        h5createGroup(filename, "bulk")
        h5write(as.vector(rownames(bulk.counts)), filename, "bulk/geneids", write.attributes = TRUE)
        h5write(as.vector(colnames(bulk.counts)), filename, "bulk/sampleids", write.attributes = TRUE)
        h5write(as.matrix(bulk.counts), filename, "bulk/data", write.attributes = TRUE)
    }
    # write bulk props if present
    if(!is.null(bulk.props)){
        h5createGroup(filename, "proportions")
        h5write(as.vector(rownames(bulk.props)), filename, "proportions/celltypeids")
        h5write(as.vector(colnames(bulk.props)), filename, "proportions/sampleids")
        h5write(as.matrix(bulk.props), filename, "proportions/data")
    }
}

read_data <- function(filename){
    require(rhdf5)
    content <- h5ls(filename, recursive = T)
    # read data that was stored using write_data
    # assume that if a group is present, all expected subgroups etc are available
    # this works as long as no external data is read (use function above)
    if("bulk" %in% content$name){
        bulk.counts <- h5read(filename, "bulk/data")
        rownames(bulk.counts) <- h5read(filename, "bulk/geneids")
        colnames(bulk.counts) <- h5read(filename, "bulk/sampleids")
    }else{
        bulk.counts <- NULL
    }
    if("proportions" %in% content$name){
        bulk.props <- h5read(filename, "proportions/data")
        rownames(bulk.props) <- h5read(filename, "proportions/celltypeids")
        colnames(bulk.props) <- h5read(filename, "proportions/sampleids")
    }else{
        bulk.props <- NULL
    }
    if("singlecell" %in% content$name){
        sc.counts <- h5read(filename, "singlecell/data")
        rownames(sc.counts) <- h5read(filename, "singlecell/geneids")
        colnames(sc.counts) <- h5read(filename, "singlecell/cellids")
        # reconstruct pheno dataframe from single vectors
        sc.pheno <- c()
        pheno.rows <- which(content$group == "/singlecell/pheno")
        pheno.names <- content[pheno.rows, "name"]
        sc.pheno <- data.frame(h5read(filename, paste("singlecell/pheno", pheno.names[1], sep = "/")))
        if(length(pheno.names) > 1){
            for(i in 2:length(pheno.names)){
                sc.pheno <- cbind(sc.pheno, h5read(filename, paste("singlecell/pheno", pheno.names[i], sep = "/")))
            }
        }
        colnames(sc.pheno) <- pheno.names
    }else{
        sc.counts <- NULL
        sc.pheno <- NULL
    }
    return(list(sc.counts = sc.counts, sc.pheno = sc.pheno, bulk.counts = bulk.counts, bulk.props = bulk.props))
}

write_misc_input <- function(genesets, algorithm.names, function.call, grouping, filename) {
    require(rhdf5)
    h5createFile(filename)
    if(!is.null(genesets)){
        h5createGroup(filename, "genesets")
        for(g in names(genesets)){
            h5write(as.vector(genesets[[g]]), filename, paste("genesets", g, sep = "/"))
        }
    }
    h5write(as.vector(algorithm.names), filename, "algorithms")
    h5write(as.vector(grouping), filename, "grouping")
    h5createGroup(filename, "function_call")
    function.args <- as.character(function.call)
    function.argnames <- names(as.list(function.call))
    h5write(as.vector(function.args), filename, "function_call/args")
    h5write(as.vector(function.argnames), filename, "function_call/argnames")
}

read_misc_input <- function(filename){
    content <- h5ls(filename)
    if("/genesets" %in% content$group){
        genesets <- list()
        for(name in content$name[which(content$group == "/genesets")]){
            genesets[[name]] <- h5read(filename, paste("genesets/",name,sep="/"))
        }
    }else{
        genesets <- NULL
    }
    algorithms <- h5read(filename, "algorithms")
    grouping <- h5read(filename, "grouping")
    function.call <- as.list(h5read(filename, "function_call/args"))
    names(function.call) <- h5read(filename, "function_call/argnames")
    function.call <- as.call(function.call)
    return(list(genesets = genesets, algorithms = algorithms, grouping = grouping, function.call = function.call))
}

write_list <- function(result.list, filename, group = NULL){
    print("calling")
    if(!file.exists(filename)) h5createFile(filename)
    if(any(sapply(result.list, function(x){is.list(x)}))){
        for(name in names(result.list)){
            if(!is.null(group)){
                groupname <- paste(group, name, sep = "/")
            }else{
                groupname <- name
            }
            h5createGroup(filename, groupname)
            if(is.list(result.list[[name]])){
                write_list(result.list[[name]], filename, groupname)
            }else{
                if(name == "bulk.props"){
                    h5write(as.matrix(result.list[[name]]), filename, paste(groupname, "data", sep = "/"))
                    h5write(as.vector(rownames(result.list[[name]])), filename, paste(groupname, "celltypeids", sep = "/"))
                    h5write(as.vector(colnames(result.list[[name]])), filename, paste(groupname, "bulkids", sep = "/"))
                }
            }
        }
    }else{
        print("down")
        print(str(result.list))
        if(!is.null(result.list[["est.props"]])){
            h5createGroup(filename, paste(group, "est.props", sep = "/"))
            h5write(as.matrix(result.list[["est.props"]]), filename, paste(group, "est.props", "data", sep = "/"))
            h5write(as.vector(rownames(result.list[["est.props"]])), filename, paste(group, "est.props", "celltypeids", sep = "/"))
            h5write(as.vector(colnames(result.list[["est.props"]])), filename, paste(group, "est.props", "bulkids", sep = "/"))
        }
        if(!is.null(result.list[["sig.matrix"]])){
            h5createGroup(filename, paste(group, "sig.matrix", sep = "/"))
            h5write(as.matrix(result.list[["sig.matrix"]]), filename, paste(group, "sig.matrix", "data", sep = "/"))
            h5write(as.vector(rownames(result.list[["sig.matrix"]])), filename, paste(group, "sig.matrix", "geneids", sep = "/"))
            h5write(as.vector(colnames(result.list[["sig.matrix"]])), filename, paste(group, "sig.matrix", "celltypeids", sep = "/"))
        }
        h5write(result.list[["name"]], filename, paste(group, "name", sep = "/"))
        h5write(result.list[["times"]], filename, paste(group, "times", sep = "/"))
    }
}

read_results <- function(filename) {
}
