#' write_data
#'
#' @param sc.counts numeric matrix, with feature as rows and scRNA-Seq
#'  profiles as columns
#' @param sc.pheno data frame, with scRNA-Seq profiles as rows, and pheno
#'  entries as rows
#' @param bulk.counts numeric matrix, with features as rows bulk measurements
#'  as columns
#' @param bulk.props numeric matrix cell type proportions, cell types as rows, 
#' bulks as columns
#' @param sub.props ???
#' @param filename string, where should the data be stored?
#'
#' @return NULL, function saves into ‘filename‘
#' another line
#' and another one
write_data <- function(sc.counts = NULL, sc.pheno = NULL, bulk.counts = NULL, bulk.props = NULL, sub.props = NULL, filename) {
  library(rhdf5)
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
  if(!is.null(sub.props)){
    h5createGroup(filename, "fine_proportions")
    h5write(as.vector(rownames(sub.props)), filename, "fine_proportions/celltypeids")
    h5write(as.vector(colnames(sub.props)), filename, "fine_proportions/sampleids")
    h5write(as.matrix(sub.props), filename, "fine_proportions/data")
  }
}