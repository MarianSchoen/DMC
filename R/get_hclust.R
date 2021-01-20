#' Cluster scRNA-Seq data hierarchically, within the input lables/cell types
#'
#' Basically, we use a 3 step approach:
#' Iteratively go through all cell types available. 
#' Withini one cell type: 
#'    1. Get the 'n.features.pre.pca' by 'feature.select.metric.FUN'
#'    2. Calculate PCA, select n PCs such that
#'     'perc.sd.in.pca' of the datas sd is explained
#'    3. Use the PCA, to embedd the data in a 2D UMAP
#'    4. Calculate the hierachical clustering in the UMAP (via R package "uwot")
#'     space
#'
#' @param sc.counts count matrix, features as rows, scRNA-Seq profiles as columns
#' @param sc.pheno data.frame. scRNA-Seq profiles as rows. Must have 
#' 'cell.type.column' and 'sample.name.column' 
#' @param cell.type.column string, column of 'sc.pheno' holding the 
#' input cell type labels. Within these entries, the clustering is done. 
#' @param sample.name.column string, column of the 'colnames(sc.counts)'
#' @param n.features.pre.pca 1 < integer , number of features selected prior to PCA
#' @param feature.select.metric.FUN function, that ranks the features of 
#' 'sc.counts'. Given a numeric vector, the function must return a single numeric.
#' We pick the top features prior to PCA. 
#' @param perc.sd.in.pca 0 < numeric <= 1. As many principal components are
#' picked, such that at least 'perc.sd.in.pca' standard deviation of the data is
#' explained
#' @param verbose logical, should information about the process be printed. 
#' @param linkage.method string, passed to hclust as 'method'. There, it says: 
#' the agglomeration method to be used. This should be (an unambiguous
#'  abbreviation of) one of "ward.D", "ward.D2", "single", "complete",
#'   "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#'    "centroid" (= UPGMC).
#'
#' @return
#'
#' @examples
get_hclust <- function(
  sc.counts,
  sc.pheno,
  cell.type.column, 
  sample.name.column, 
  n.features.pre.pca = 4e3,
  feature.select.metric.FUN = sd,
  perc.sd.in.pca = 0.8, 
  verbose = TRUE, 
  linkage.method = "average"
) {
  
  #library(uwot)
  
  # few basic safety checks (all default values are just ideas):
  if (!is.numeric(n.features.pre.pca)) {
    n.features.pre.pca <- round(nrow(sc.counts) * 0.1)
  } else {
    if (n.features.pre.pca > nrow(sc.counts)) {
      n.features.pre.pca <- nrow(sc.counts)
    }
    if (n.features.pre.pca <= 1) {
      n.features.pre.pca <- nrow(sc.counts)
    }
  }
  
  if (!is.numeric(perc.sd.in.pca)) {
    perc.sd.in.pca <- 0.8
  }
  
  if (!is.function(feature.select.metric.FUN)) {
    feature.select.metric.FUN <- sd
  } else {
    # test whether feature.select.metric.FUN returns a
    # single value given a numeric vector
    tmp <- feature.select.metric.FUN(1:10)
    if (length(tmp) != 1 || !is.numeric(tmp)) {
      feature.select.metric.FUN <- sd
    }
  }
  
  if ( !cell.type.column %in% colnames(sc.pheno) ) {
    stop("In get_hclust: 'cell.type.column' not in 'colnames(sc.pheno)'")
  }
  
  if ( !sample.name.column %in% colnames(sc.pheno) ) {
    stop("In get_hclust: 'sample.name.column' not in 'colnames(sc.pheno)'")
  }
  
  
  cell.types <- unique(as.character(sc.pheno[[cell.type.column]]))
  # select features using the provided function:
  if(verbose){
    print(
      paste0(
        "Detected ", length(cell.types), " cell.types"
      )
    )
    print(
      paste0(
        "Starting to iteratively cluster them at " , Sys.time()
        )
      )
  }
  
  clust.obj <- list()
  
  for( cell.type in cell.types){
    if( verbose) {
        print(paste0("Starting with ", cell.type, " at ", Sys.time()))
      }
    profiles.of.this.type <- which(
      as.character(sc.pheno[[cell.type.column]]) == cell.type
      )
    
    profile.names.of.this.type <- sc.pheno[[sample.name.column]][profiles.of.this.type]
    
    intersected.profile.names <- intersect(
      profile.names.of.this.type, 
      colnames(sc.counts)
      )
    
    if( verbose ) {
      print(
        paste0(
          "There are ", length(intersected.profile.names), " profiles"
          )
        )
    }    
    
    if(length(intersected.profile.names) == 0){
      warning(
        paste0(
          "In 'get_hclust': can't find columns of type ", cell.type, 
          " in 'sc.counts' if intersecting ",
          "sc.pheno[[sample.name.column]]' and 'colnames(sc.counts)'"
          )
      )
      clust.obj[[cell.type]] <- NA
      next
    }
    
    tmp.counts <- sc.counts[, profile.names.of.this.type]
    
    sd.per.gene <- apply(tmp.counts, 1, feature.select.metric.FUN)
    names(sd.per.gene) <- rownames(tmp.counts)
    sorted.sds <- sort(sd.per.gene, decreasing = TRUE)
    preselected.features <- names(sorted.sds)[1:n.features.pre.pca]
    
    tmp.pca <- prcomp(
      x = t(Matrix::as.matrix(tmp.counts[which(rownames(tmp.counts) %in% preselected.features), ]))
      # , rank. = 100 # this parameter is useless, it calculates all pcs, and slices the output matrix
      , scale = FALSE
      , center = FALSE
    )
    
    perc.sd <- tmp.pca$sdev/sum(tmp.pca$sdev)
    
    cumulative.perc.sd <- sapply(
      X = 1:ncol(tmp.counts)
      , FUN = function(pos){
        sum(perc.sd[1:pos])
      })
    n.components <- which(cumulative.perc.sd > perc.sd.in.pca)[1]
    
    tmp.umaps <- uwot::umap(
      X = tmp.pca$x[, 1:n.components]
      , n_neighbors = ifelse(
        test = length(profiles.of.this.type) > 60 # an arbitrary number
        , yes = 60
        , no = length(profiles.of.this.type)
      )
      , scale = FALSE
      , pca = NULL # according to the documentation, this means "no pca"
    )
    
    clust.obj[[cell.type]] <- hclust(
      d = dist(tmp.umaps)
      , method = linkage.method
    )
  }
  return(clust.obj)
}
