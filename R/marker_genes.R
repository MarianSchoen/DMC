#' select marker genes according to Baron et al., 2016
#' 
#' @param exprs matrix containing single cell profiles as columns
#' @param pheno phenotype data corresponding to the expression matrix.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param sig.types character vector containing cell types for which 
#' marker genes should be selected. Default is all of them.
#' @return list containing a vector of marker genes for each cell type
#' @example marker_genes(sc.exprs, sc.pheno)

marker_genes <- function(exprs, pheno, sig.types = NULL){
  # error checking
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  if (!all(sig.types %in% pheno[, "cell_type"])) {
    stop("Not all cell types in 'sig.types' are present in pheno data.")
  }

  exprs <- scale_to_count(exprs)
  if(is.null(sig.types)) sig.types <- unique(pheno[,"cell_type"])
  
  # 1) not sure this interpretation is right, but I cannot think of something more reasonable
  geneSums <- rowSums(exprs)
  
  # calculate sum per gene per cell type and the fano factor (actually just select by variance) per gene per cell type
  sums.per.celltype <- matrix(0, nrow = nrow(exprs), ncol = length(sig.types))
  colnames(sums.per.celltype) <- sig.types
  rownames(sums.per.celltype) <- rownames(exprs)
  fano.per.celltype <- matrix(0, nrow = nrow(exprs), ncol = length(sig.types))
  colnames(fano.per.celltype) <- sig.types
  rownames(fano.per.celltype) <- rownames(exprs)
  for(t in sig.types){
    sums.per.celltype[,t] <- apply(exprs[,which(pheno[,"cell_type"] == t),drop=F], 1, sum) / geneSums
    fano.per.celltype[,t] <- apply(exprs[,which(pheno[,"cell_type"] == t),drop=F], 1, var)
  }
  
  # take only genes with variance above the median
  fano.per.celltype[is.na(fano.per.celltype)] <- 0
  
  # only highly expressed genes
  sums.above.threshold <- apply(sums.per.celltype, 2, function(x) x>0.5)
  
  # genes with high variance
  factors.above.theshold <- apply(fano.per.celltype,2,function(x){x > summary(x)[3]})
  
  # for each cell type, take the intersect of the first criteria
  genes.per.celltype <- list()
  for(t in sig.types){
    temp.genes <- rownames(exprs)[intersect(which(sums.above.threshold[,t]), which(factors.above.theshold[,t]))]
    if(length(temp.genes) > 0){
	    genes.per.celltype[[t]] <- temp.genes
    }else{
	    genes.per.celltype[t] <- list(NULL)
  	}
  }
  
  # do a ks test for each of these genes for each group against all others combined?
  for(t in sig.types){
    grouping <- ifelse(pheno[,"cell_type"] == t, TRUE, FALSE)
    genes <- genes.per.celltype[[t]]
    if(length(genes)>0){
      p.vals <- apply(exprs[genes,,drop=FALSE], 1, function(x){ks.test(x[grouping], x[!grouping], alternative = "two.sided")$p.value})
      if(length(genes[p.vals < 10^(-5)]) > 0) genes.per.celltype[[t]] <- genes[p.vals < 10^(-5)]
    }
  }
  #print(str(genes.per.celltype))
  return(genes.per.celltype)
}
