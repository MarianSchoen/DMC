#' add subtype pheno 
#'
#' @param sc.counts 
#' @param sc.pheno 
#' @param cell.type.column 
#' @param sample.name.column 
#' @param new.subtype.column 
#' @param hclust.obj 
#' @param avg.profiles.per.subcluster.vec 
#' @param verbose 
#' @param ... 
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