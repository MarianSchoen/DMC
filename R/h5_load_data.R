#' read_data
#'
#' @param filename string, which file should be read in?
#'
#' @return list with: \cr
#'  - sc.counts: numeric matrix, features as rows, scRNA-Seq profiles as colums\cr
#'  - sc.pheno: dataframe with scRNA-Seq profiles as rows, and multiple pheno
#'   columns \cr
#'  - bulk.counts: numeric matrix, features as rows, bulk measurements as 
#'  columns\cr
#'  - bulk.props: numeric matrix, cell types as columns, bulks as columns
read_data <- function(filename){
  library(rhdf5)
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
  if("fine_proportions" %in% content$name){
    sub.props <- h5read(filename, "fine_proportions/data")
    rownames(sub.props) <- h5read(filename, "fine_proportions/celltypeids")
    colnames(sub.props) <- h5read(filename, "fine_proportions/sampleids")
  }else{
    sub.props <- NULL
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
    if(any(pheno.names == "sample.name")){
	    rownames(sc.pheno) <- sc.pheno[["sample.name"]]
    }
  }else{
    sc.counts <- NULL
    sc.pheno <- NULL
  }
  return(list(sc.counts = sc.counts, sc.pheno = sc.pheno, bulk.counts = bulk.counts, bulk.props = bulk.props, sub.props = sub.props))
}
