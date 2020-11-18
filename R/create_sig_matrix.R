#' create a signature matrix for usage with CIBERSORT or DeconRNASeq
#' according to description in Newman et. al
#'
#' @param exprs non-negative numeric matrix, 
#' @param pheno data.frame, phenodata corresponding to
#' expression data, same ordering of samples/cells; cell type
#' label column must be named `cell.type.column`
#' @param exclude.celltypes character vector, defaults to
#' c("malignant", "not_annotated", "unassigned"). Cells
#' with these labels are not included in the signature matrix
#' @param max.genes integer, maximum number of genes to be selected
#' for each cell type; the actual number of genes in the
#' signature matrix will be larger. If no value is given,
#' max.genes = NULL and all genes will be selected
#' @param optimize boolean, default TRUE. If TRUE, the
#' optimal number of genes to be selected for each cell type will be
#' determined between 1 and 'max.genes' by optimizing the
#' condition number of the signature matrix.
#' If FALSE, for each cell type the top 'max.genes' genes will be taken.
#' @param split.data logical, should the training data be split for
#' reference profile creation and optimization? default: TRUE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? 
#' @return numeric matrix, holding reference profiles in its column, 
#' features in its rows

create_sig_matrix <- function(
  exprs,
  pheno,
  exclude.celltypes = NULL,
  max.genes = NULL,
  optimize = TRUE,
  split.data = TRUE,
  cell.type.column = "cell_type"
  ) {
  # parameter checks
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  if (!is.null(max.genes) && max.genes == 0) {
        max.genes <- NULL
  }
  rownames(pheno) <- colnames(exprs)

  # exclude specified cell types
  if (!is.null(exclude.celltypes)) {
    to.exclude <- which(pheno[, cell.type.column] %in% exclude.celltypes)
    if(length(to.exclude) > 0){
      exprs <- exprs[, -to.exclude, drop = F]
      pheno <- pheno[-to.exclude, , drop = F]
    }
  }

  # make sure that not more genes than available are selected
  if (is.null(max.genes)) {
    max.genes <- floor(nrow(exprs) / length(unique(pheno[, cell.type.column])))
  }

  
  # only split data if there are at least 3 samples for each cell type
  # otherwise the testing in the sig. matrix building will fail
  type.counts <- table(pheno[, cell.type.column])

  if (split.data && !any(type.counts < 3)) {
    # create reference profiles from 30% of the cells of each type
    cell.types <- pheno[, cell.type.column]
    names(cell.types) <- colnames(exprs)

    sample.X <- DTD::sample_random_X(
      included.in.X = unique(cell.types),
      pheno = cell.types,
      expr.data = Matrix::as.matrix(exprs),
      percentage.of.all.cells = 0.3,
      normalize.to.count = TRUE
    )
    sig.matrix <- sample.X$X.matrix
    samples.to.remove <- sample.X$samples.to.remove
    rm(sample.X)
    exprs <- exprs[, -which(colnames(exprs) %in% samples.to.remove), drop = F]
    pheno <- pheno[-which(rownames(pheno) %in% samples.to.remove), , drop = F]
  }

  # for each cell type test against all others for DEG
  # using two-sided t-test
  deg.per.type <- list()
  for (t in unique(pheno[, cell.type.column])) {
    labs <- ifelse(pheno[, cell.type.column] == t, 0, 1)
    no.var.genes <- unique(
      c(
        which(apply(exprs[, which(labs == 0), drop = F], 1, var) == 0),
        which(apply(exprs[, which(labs == 1), drop = F], 1, var) == 0)   
      )
    )
    if (length(no.var.genes) > 0) {
      temp.exprs <- exprs[-no.var.genes, , drop = F]
      t.test.result <- multtest::mt.teststat(
        X = temp.exprs,
        classlabel = labs,
        test = "t"
      )

      # calc p-values and adjust for multiple testing, in accordance with Newman et al. 2015 significant if q < 0.3
      p.vals <- 2 * pt(abs(t.test.result), length(labs) - 2, lower.tail = FALSE)
      q.vals <- p.adjust(p.vals, "BH")
      sig.entries <- which(q.vals < 0.3)
      sig.genes <- rownames(temp.exprs)[sig.entries]
      # catch possible errors related to sig.genes
      if(any(is.na(sig.genes))){
        sig.entries <- sig.entries[!is.na(sig.genes)]
        sig.genes <- sig.genes[!is.na(sig.genes)]
      }
      if(!length(sig.genes) > 0) break
      temp.exprs <- temp.exprs[sig.entries, ]
      
      # the genes are ordered by decreasing fold changes compared to other cell subsets
      fold.changes <- 
        log2(Matrix::rowMeans(
          temp.exprs[, which(labs == 0), drop = F]
        )) - 
        log2(Matrix::rowMeans(
          temp.exprs[, which(labs == 1), drop = F]
        ))
    } else {
      t.test.result <- multtest::mt.teststat(
        X = exprs,
        classlabel = labs,
        test = "t"
      )
      p.vals <- 2 * pt(abs(t.test.result), length(labs) - 2, lower.tail = FALSE)
      q.vals <- p.adjust(p.vals, "BH")
      sig.entries <- which(q.vals < 0.3)
      sig.genes <- rownames(exprs)[sig.entries]
      
      # catch possible errors related to sig.genes
      if(any(is.na(sig.genes))){
        sig.entries <- sig.entries[!is.na(sig.genes)]
        sig.genes <- sig.genes[!is.na(sig.genes)]
      }
      if(!length(sig.genes) > 0) break
      
      # calculate fold changes
      fold.changes <- 
        log2(Matrix::rowMeans(
          exprs[sig.entries, which(labs == 0), drop = F]
        )) - 
        log2(Matrix::rowMeans(
          exprs[sig.entries, which(labs == 1), drop = F]
        ))
    }
    # add genes for each type ordered by decreasing fold change
    deg.per.type[[t]] <- sig.genes[order(abs(fold.changes), decreasing = TRUE)]
  }

  # reduce to one reference profile per cell type
  ref.profiles <- matrix(
    0,
    nrow = nrow(exprs),
    ncol = length(unique(pheno[, cell.type.column]))
  )
  colnames(ref.profiles) <- unique(pheno[, cell.type.column])
  rownames(ref.profiles) <- rownames(exprs)
  
  for (t in colnames(ref.profiles)) {
    ref.profiles[, t] <- Matrix::rowMeans(
      exprs[, which(pheno[, cell.type.column] == t), drop = F]
    )
  }
  
  # again following Newman et al.:
  # take top g genes for every cell type, create signature matrices
  # choose the gene set that minimizes condition number
  limit <- min(
    max(sapply(deg.per.type, length)), max.genes
  )

  if (optimize) {
    cond.nums <- rep(-1, times = limit)
    for (g in 1:limit) {
      all.genes <- unique(unlist(
        sapply(deg.per.type, function(sub.genes, lim = g){
          sub.genes[1:min(length(sub.genes), lim)]
        })
      ))

      # mt.multtest might return NA
      if (any(is.na(all.genes))) {
        all.genes <- all.genes[-which(is.na(all.genes))]
      }

      # estimate condition number of reduced matrix
      cond.nums[g] <- kappa(ref.profiles[all.genes, ,drop=F], exact = FALSE)
    }
    optimal.g <- which.min(cond.nums)
  } else {
    optimal.g <- limit
  }

  # create gene list with the optimal g
  optimal.genes <- unique(unlist(
    sapply(deg.per.type, function(sub.genes, opt.g = optimal.g) {
      sub.genes[1:min(opt.g, length(sub.genes))]
    })
  ))

  # deal with possible NAs again
  if (any(is.na(optimal.genes))) {
    optimal.genes <- optimal.genes[-which(is.na(optimal.genes))]
  }

  # once again depending on whether split.data is true and possible
  if (split.data && !any(type.counts < 3)) {
    ref.mat <- sig.matrix[optimal.genes, ]
    rm(sig.matrix)
  } else {
    ref.mat <- ref.profiles[optimal.genes, ]
    rm(ref.profiles)
  }

  # remove any duplicate genes that might be in the matrix
  if(any(duplicated(rownames(ref.mat)))){
    warning("Found duplicates in reference matrix")
    ref.mat <- ref.mat[-which(duplicated(rownames(ref.mat))),]
  }
  return(ref.mat)
}
