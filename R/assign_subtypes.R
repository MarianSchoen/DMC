#' assign_subtypes
#'
#' @param sc.counts non-negative numeric matrix with features as rows, and 
#' scRNA-Seq profiles as columns. 'ncol(sc.counts)' must equal 'nrow(sc.pheno)'
#' @param sc.pheno data frame with scRNA-Seq profiles as rows, and pheno entries
#' in columns. 'nrow(sc.pheno)' must equal 'ncol(sc.counts)' 
#' @param sub.list #TODO I don't understand
#' @param celltypecol string, column name in 'sc.pheno' that holds the cell type
#' per scRNA-Seq profile
#' @param ... additional parameters that get passed to Rtsne()
#'
#' @return
#' @export
#'
#' @examples
assign_subtypes <- function(
  sc.counts
  , sc.pheno
  , sub.list
  , celltypecol = "cell_type"
  , ...
  ){
  require(Rtsne)
  # parameter checks
  if(!celltypecol %in% names(sc.pheno)) stop("celltype column not in data frame")
  if("subtype" %in% names(sc.pheno)) stop("subtype column already present")
  if(ncol(sc.counts) != nrow(sc.pheno)) stop("expression and pheno data do not match")
  tsne.embed <- Rtsne(t(sc.counts), ...)
  # add default subtype 1
  sc.pheno <- cbind(sc.pheno, subtype = rep(1, nrow(sc.pheno)))
  for(ct in names(sub.list)){
    if(!ct %in% sc.pheno[[celltypecol]]) {
      warning(paste(ct, ": No such cell type"))
      next
    }
    ct.indices <- which(sc.pheno[[celltypecol]] == ct)
    km.clust <- kmeans(x = tsne.embed$Y[ct.indices, ], min(sub.list[[ct]], length(ct.indices)))
    sc.pheno[ct.indices, "subtype"] <- km.clust$cluster
  }
  return(list(sc.pheno = sc.pheno, tsne.embed = tsne.embed))
}