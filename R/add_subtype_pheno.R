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
#' @param sc.counts count matrix, features as rows, scRNA-Seq profiles as columns
#' @param sc.pheno data.frame. scRNA-Seq profiles as rows. Must have 
#' 'cell.type.column' and 'sample.name.column' 
#' @param cell.type.column string, column of 'sc.pheno' holding the 
#' input cell type labels. Within these entries, the clustering is done. 
#' @param sample.name.column string, column of the 'colnames(sc.counts)'
#' @param new.subtype.column string, pattern of the new column, where the subtype
#' information is stored
#' @param hclust.obj 'hclust' object, or NA. If you already calculated a hclust, 
#' object pass it here. Alternatively, we call the get_hclust function. 
#' @param avg.profiles.per.subcluster.vec numeric vector, indicates how many new 
#' subtype pheno entries are generated 
#' (=> 'length(avg.profiles.per.subcluster.vec)') and - in average - how many
#' profiles are in a subtype.
#' @param verbose logical, should information about the process be printed. 
#' @param ... arguments that are passed to get_hclust
#'
#' @return
#'
#' @examples
add_subtype_pheno <- function(
  sc.counts = cll.exprs,
  sc.pheno = cll.pheno,
  cell.type.column = "cell_type", 
  sample.name.column = "sample.name", 
  new.subtype.column = "subtype",
  hclust.obj = NA,
  avg.profiles.per.subcluster.vec = seq(from = 100, to = 1.1e3, by = 500),
  verbose = TRUE,
  ...
){
  if(! is.list(hclust.obj) ){
    hclust.obj <- get_hclust(
      sc.counts = sc.counts
      , sc.pheno = sc.pheno
      , cell.type.column = cell.type.column
      , sample.name.column = sample.name.column
      , verbose = verbose
      , ...
    )
  }
  for(avg.profiles.per.subcluster in avg.profiles.per.subcluster.vec){
    new.entry.name <- paste0(
      new.subtype.column, ".avg.", avg.profiles.per.subcluster
    )
    sc.pheno[[new.entry.name]] <- NA
    
    for(cell.type in unique(sc.pheno[[cell.type.column]])){
      profiles.pos <- which(sc.pheno[[cell.type.column]] == cell.type)
      
      n.new.clusters <- ceiling(length(profiles.pos)/avg.profiles.per.subcluster)
      if(n.new.clusters < 2){
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
  
  return(sc.pheno)
}
