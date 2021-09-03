#' add subtype pheno
#'
#' Use a hierarchical clustering of the input labels,
#' to generate finer subtype labels.
#' The basic steps are:
#'    1. use the get_hclust function to hierarchically cluster the data
#'     within the cell types. (or use the provided hclust.obj)
#'    2. cluster the data multiple times, with different depth
#'    (as indicated by 'avg.profiles.per.subcluster.vec')
#'    3. provide the clustering labels as columns in sc.pheno, and return it.
#'
#' @param sc.counts count matrix, features as rows,
#' scRNA-Seq profiles as columns
#' @param sc.pheno data.frame. scRNA-Seq profiles as rows. Must have
#' 'cell.type.column' and 'sample.name.column'
#' @param cell.type.column string, column of 'sc.pheno' holding the
#' input cell type labels. Within these entries, the clustering is done.
#' @param sample.name.column string, column of the 'colnames(sc.counts)'
#' @param new.subtype.column string, pattern of the new column,
#' where the subtype information is stored
#' @param hclust.obj 'hclust' object, or NA. If you already calculated a hclust,
#' object pass it here. Alternatively, we call the get_hclust function.
#' @param n.clusters integer vector of clustering depths (number of subclusters
#' created for each cell type), default c(1, 2, 4, 8).
#' This means that in the finest clustering, each celltype will be split in
#' 8 subtypes, in the next step each will be split in 4 subtypes, ...
#' @param verbose logical, should information about the process be printed.
#' @param ... arguments that are passed to get_hclust
#'
#' @return data.frame sc.pheno, with added subtype column
#'
#' @examples
add_subtype_pheno <- function(
  sc.counts,
  sc.pheno,
  cell.type.column = "cell_type",
  sample.name.column = "sample.name",
  new.subtype.column = "subtype",
  hclust.obj = NA,
  n.clusters = c(1, 2, 4, 8),
  verbose = TRUE,
  ...
) {
	cat("Adding Subtypes\n")
  if (! is.list(hclust.obj)) {
    hclust.obj <- get_hclust(
      sc.counts = sc.counts,
      sc.pheno = sc.pheno,
      cell.type.column = cell.type.column,
      sample.name.column = sample.name.column,
      verbose = verbose,
      ...
    )
  }

  for (n.new.clusters in n.clusters) {
    new.entry.name <- paste0(
      new.subtype.column, ".frac.", n.new.clusters
    )
    sc.pheno[[new.entry.name]] <- NA

    for (cell.type in unique(sc.pheno[[cell.type.column]])) {
      profiles.pos <- which(sc.pheno[[cell.type.column]] == cell.type)

      if (n.new.clusters < 2) {
        sc.pheno[[new.entry.name]][profiles.pos] <- cell.type
        next
      }

      subclusters <- cutree(
        tree = hclust.obj[[cell.type]]
        , k = n.new.clusters
      )

      sc.pheno[[new.entry.name]][profiles.pos] <- paste0(
        cell.type, ".", subclusters
      )
    }
  }

  cat("Done\n")
  rm(hclust.obj)
  gc()
  return(sc.pheno)
}
