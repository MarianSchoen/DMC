#' simulate subtypes of given cell types using PCA and hierarchical clustering
#'
#' @param sc.counts non-negative numeric matrix with features as rows, and 
#' scRNA-Seq profiles as columns. 'ncol(sc.counts)' must equal 'nrow(sc.pheno)'
#' @param sc.pheno data frame with scRNA-Seq profiles as rows, and pheno entries
#' in columns. 'nrow(sc.pheno)' must equal 'ncol(sc.counts)' 
#' @param sub.list named list containing the number of subtypes to be simulated for each cell type
#' @param celltypecol string, column name in 'sc.pheno' that holds the cell type
#' per scRNA-Seq profile
#'
#' @return list containing two entries:  
#' 
#' 1) sc.pheno - input pheno data with added column 'subtype'  
#' 
#' 2) tsne.embed - t-SNE embedding (output of Rtsne) that was used to generate
#' the subtypes
#'
#' @examples
assign_subtypes <- function(
  sc.counts
  , sc.pheno
  , sub.list
  , celltypecol = "cell_type"
  ){
  # parameter checks
  if(!celltypecol %in% names(sc.pheno)) stop("celltype column not in data frame")
  if("subtype" %in% names(sc.pheno)) {
    warning("subtype column already present. returning data frame unchanged")
    return(list(sc.pheno = sc.pheno, tsne.embed = NULL))
  }
  if(ncol(sc.counts) != nrow(sc.pheno)) stop("number of columns in sc.counts and rows in sc.pheno do not match")
  if(!is.list(sub.list)){
    stop("sub.list must be a list")
  }

  # calculate PCA
  counts.pca <- prcomp(sc.counts)$x[, 1:50]
  # hierarchical clustering
  dists <- as.matrix(dist(counts.pca, "euclidean"))

  # add default subtype 1
  sc.pheno <- cbind(sc.pheno, subtype = rep(1, nrow(sc.pheno)))

  # if tsne worked, cluster each celltype in tsne-space and 
  # assign each cluster a subtype
  if(!class(tsne.embed) == "try-error"){
    for(ct in names(sub.list)){
      if(!ct %in% sc.pheno[[celltypecol]] || sub.list[[ct]] < 2) {
        warning(paste(ct, ": No such cell type"))
        next
      }
      ct.indices <- which(sc.pheno[[celltypecol]] == ct)

      # hierarchical clustering on pca
      temp.dists <- as.dist(dists[ct.indices, ct.indices])
      tree <- hclust(temp.dists, "average")
      clustering <- cutree(tree, k = sub.list[[ct]])

      sc.pheno[ct.indices, "subtype"] <- clustering
    }
  }
  
  return(list(sc.pheno = sc.pheno))
}