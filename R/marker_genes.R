#' select marker genes following description in Baron et al., 2016
#' 
#' @param exprs matrix containing single cell profiles as columns
#' @param pheno phenotype data corresponding to the expression matrix.
#' Has to contain single cell labels in a column named `cell.type.column`
#' @param sig.types character vector containing cell types for which 
#' marker genes should be selected. Default is all of them.
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information?
#' @return list containing a vector of marker genes for each cell type
#' @example marker_genes(sc.exprs, sc.pheno)

marker_genes <- function(exprs, pheno, sig.types = NULL, cell.type.column = "cell_type"){
  # parameter checks
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  if (!all(sig.types %in% pheno[, cell.type.column])) {
    stop("Not all cell types in 'sig.types' are present in pheno data.")
  }
  if(is.null(sig.types)) sig.types <- unique(pheno[,cell.type.column])
  rownames(pheno) <- colnames(exprs)

  # scale to total count number for testing
  exprs <- scale_to_count(exprs)

  # calculate sum per gene per cell type and the fano 
  # factor (select by variance) per gene per cell type
  # see Baron et al., 2016 (Methods)
  geneSums <- Matrix::rowSums(exprs)

  # create matrices for variance and sum per celltype and gene
  sums.per.celltype <- matrix(0, nrow = nrow(exprs), ncol = length(sig.types))
  colnames(sums.per.celltype) <- sig.types
  rownames(sums.per.celltype) <- rownames(exprs)
  fano.per.celltype <- matrix(0, nrow = nrow(exprs), ncol = length(sig.types))
  colnames(fano.per.celltype) <- sig.types
  rownames(fano.per.celltype) <- rownames(exprs)

  for(t in sig.types){
    sums.per.celltype[,t] <- Matrix::rowSums(
      exprs[,which(pheno[,cell.type.column] == t),drop=F]
    ) / geneSums
    fano.per.celltype[,t] <- apply(
      exprs[,which(pheno[,cell.type.column] == t),drop=F], 
      1, 
      function(x) {var(x) / mean(x)})
  }
  
  fano.per.celltype[is.na(fano.per.celltype)] <- 0
  
  # only highly expressed genes
  sums.above.threshold <- apply(sums.per.celltype, 2, function(x) x>0.5)
  
  # genes with high variance
  factors.above.theshold <- apply(fano.per.celltype,2,function(x){x > summary(x)[3]})
  
  # for each cell type, take the intersect of the first criteria
  genes.per.celltype <- list()
  for(t in sig.types){
    temp.genes <- rownames(exprs)[intersect(
        which(sums.above.threshold[,t]), 
        which(factors.above.theshold[,t])
      )]
    if(length(temp.genes) > 0){
	    genes.per.celltype[[t]] <- temp.genes
    }else{
	    genes.per.celltype[t] <- list(NULL)
  	}
  }
  
  # calculate p-values for the genes that passed the first two criteria 
  # (using ks test) and keep only those with p<10^-5 (Baron et al. 2016)
  for(t in sig.types){
    grouping <- ifelse(pheno[, cell.type.column] == t, TRUE, FALSE)
    genes <- genes.per.celltype[[t]]
    if(length(genes)>0){
      p.vals <- apply(
        exprs[genes,,drop=FALSE], 
        1, 
        function(x){
          ks.test(x[grouping], x[!grouping], alternative = "two.sided")$p.value
        }
      )
      top.genes <- genes[p.vals < 10^(-5)]
      if(length(top.genes) > 0){
        genes.per.celltype[[t]] <- top.genes
      } 
    }
  }

  return(genes.per.celltype)
}
