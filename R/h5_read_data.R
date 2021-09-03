#' read data conforming to the format used by \link{benchmark} for data storage
#'
#' @param filename string, which file should be read in?
#' @return list with: \cr
#'  - sc.counts: numeric matrix, features as rows,
#'    scRNA-Seq profiles as colums\cr
#'  - sc.pheno: dataframe with scRNA-Seq profiles as rows, and multiple pheno
#'   columns \cr
#'  - bulk.counts: numeric matrix, features as rows, bulk measurements as
#'  columns\cr
#'  - bulk.props: numeric matrix, cell types as rows, bulks as columns\cr
#'  - bulk.pheno: dataframe containing pheno data for bulks in columns,
#'  bulks as rows
#'@export

read_data <- function(filename) {
  # parameter check
  if (!file.exists(filename)) {
    stop(paste("Could not find file ", filename, sep = ""))
  }

  content <- rhdf5::h5ls(filename, recursive = TRUE)
  # read data that was stored using write_data
  # assume that if a group is present, all expected subgroups etc are available
  # this works as long as no external data is read (use function above)
  if ("bulk" %in% content$name) {
    bulk.counts <- Matrix(rhdf5::h5read(filename, "bulk/data"), sparse = TRUE)
    class(bulk.counts) <- "dgCMatrix"
    geneids <- rhdf5::h5read(filename, "bulk/geneids")
    sampleids <- rhdf5::h5read(filename, "bulk/sampleids")
    if (nrow(bulk.counts) == length(geneids) &&
        ncol(bulk.counts) == length(sampleids)) {
	    rownames(bulk.counts) <- geneids
	    colnames(bulk.counts) <- sampleids
    }else{
	    colnames(bulk.counts) <- geneids
	    rownames(bulk.counts) <- sampleids
	    bulk.counts <- t(bulk.counts)
    }
  }else{
    bulk.counts <- NULL
  }
  if ("/bulk/pheno" %in% content$group) {
    # reconstruct pheno dataframe from single vectors
    bulk.pheno <- c()
    pheno.rows <- which(content$group == "/bulk/pheno")
    pheno.names <- content[pheno.rows, "name"]
    bulk.pheno <- data.frame(
      rhdf5::h5read(filename, paste("bulk/pheno", pheno.names[1], sep = "/"))
    )
    if (length(pheno.names) > 1) {
      for (i in 2:length(pheno.names)) {
        bulk.pheno <- cbind(
          bulk.pheno,
          rhdf5::h5read(
            filename, paste("bulk/pheno", pheno.names[i], sep = "/")
          )
        )
      }
    }
    colnames(bulk.pheno) <- pheno.names
  }else{
    bulk.pheno <- NULL
  }
  if ("proportions" %in% content$name) {
    bulk.props <- rhdf5::h5read(filename, "proportions/data")
    celltypeids <- rhdf5::h5read(filename, "proportions/celltypeids")
    sampleids <- rhdf5::h5read(filename, "proportions/sampleids")

    # make sure the matrix dimensions are in the right order
    if (length(sampleids) == ncol(bulk.props) &&
       length(celltypeids) == nrow(bulk.props)) {
	    rownames(bulk.props) <- celltypeids
	    colnames(bulk.props) <- sampleids
    }else{
	    rownames(bulk.props) <- sampleids
	    colnames(bulk.props) <- celltypeids
	    bulk.props <- t(bulk.props)
    }
  }else{
    bulk.props <- NULL
  }

  if ("singlecell" %in% content$name) {
    sc.counts <- Matrix(
      rhdf5::h5read(filename, "singlecell/data"),
      sparse = TRUE
    )
    class(sc.counts) <- "dgCMatrix"
    geneids <- rhdf5::h5read(filename, "singlecell/geneids")
    cellids <- rhdf5::h5read(filename, "singlecell/cellids")

    # make sure the matrix dimensions are in the right order
    if (nrow(sc.counts) == length(geneids) &&
       ncol(sc.counts) == length(cellids)) {
	    rownames(sc.counts) <- geneids
	    colnames(sc.counts) <- cellids
    }else{
	    colnames(sc.counts) <- geneids
	    rownames(sc.counts) <- cellids
	    sc.counts <- t(sc.counts)
    }
  }else{
    sc.counts <- NULL
  }
  if ("/singlecell/pheno" %in% content$group) {
    # reconstruct pheno dataframe from single vectors
    sc.pheno <- c()
    pheno.rows <- which(content$group == "/singlecell/pheno")
    pheno.names <- content[pheno.rows, "name"]
    sc.pheno <- data.frame(
      rhdf5::h5read(
        filename, paste("singlecell/pheno", pheno.names[1], sep = "/")
      )
    )
    if (length(pheno.names) > 1) {
      for (i in 2:length(pheno.names)) {
        sc.pheno <- cbind(
          sc.pheno,
          rhdf5::h5read(
            filename, paste("singlecell/pheno", pheno.names[i], sep = "/")
          )
        )
      }
    }
    colnames(sc.pheno) <- pheno.names
  }else{
    sc.pheno <- NULL
  }
  return(list(
    sc.counts = sc.counts,
    sc.pheno = sc.pheno,
    bulk.counts = bulk.counts,
    bulk.props = bulk.props,
    bulk.pheno = bulk.pheno
  ))
}
