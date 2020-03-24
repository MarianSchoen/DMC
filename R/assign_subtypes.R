assign_subtypes <- function(sc.exprs, sc.pheno, sub.list, celltypecol = "cell_type", ...){
  require(Rtsne)
  # parameter checks
  if(!celltypecol %in% names(sc.pheno)) stop("celltype column not in data frame")
  if("subtype" %in% names(sc.pheno)) stop("subtype column already present")
  if(ncol(sc.exprs) != nrow(sc.pheno)) stop("expression and pheno data do not match")
  tsne.embed <- Rtsne(t(sc.exprs), ...)
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