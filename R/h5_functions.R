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