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

  # exclude genes with zero variance
  if (any(apply(exprs, 1, var) == 0)) {
    exprs <- exprs[-which(apply(exprs, 1, var) == 0), , drop = F]
  }

  # make sure that not more genes than available are selected
  if (is.null(max.genes)) {
    max.genes <- floor(nrow(exprs) / length(unique(pheno[, cell.type.column])))
  }

  # only split data if there are at least 3 samples for each cell type
  # otherwise the testing in the sig. matrix building will fail
  type.counts <- c()
  for (t in unique(pheno[, cell.type.column])) {
    type.counts <- c(type.counts, length(which(pheno[, cell.type.column] == t)))
  }

  if (split.data && !any(type.counts < 3)) {
    # create reference profiles from 10% of the cells of each type
    cell.types <- pheno[, cell.type.column]
    names(cell.types) <- colnames(exprs)

    sample.X <- DTD::sample_random_X(
      included.in.X = unique(cell.types),
      pheno = cell.types,
      expr.data = exprs,
      percentage.of.all.cells = 0.3,
      normalize.to.count = TRUE
    )
    sig.matrix <- sample.X$X.matrix
    samples.to.remove <- sample.X$samples.to.remove
    exprs <- exprs[, -which(colnames(exprs) %in% samples.to.remove), drop = F]
    pheno <- pheno[-which(rownames(pheno) %in% samples.to.remove), , drop = F]
  }

  # for each cell type test against all others for DEG
  # using two-sided t-test
  deg.per.type <- list()
  for (t in unique(pheno[, cell.type.column])) {
    labs <- ifelse(pheno[, cell.type.column] == t, 0, 1)
    no.var.genes <- which(
      apply(exprs[, which(labs == 0), drop = F], 1, var) == 0
    )
    if (length(no.var.genes) > 0) {
      t.test.result <- multtest::mt.teststat(
        X = exprs[-no.var.genes, , drop = F],
        classlabel = labs,
        test = "t"
      )

      # calc p-values and adjust for multiple testing, in accordance with Newman et al. 2015 significant if q < 0.3
      p.vals <- 2 * pt(abs(t.test.result), length(labs) - 2, lower.tail = FALSE)

      # apparently an error occurs when qvalue calls pi0est, which calls smooth.spline; input does not contain
      # missing or infinite values, thereforce catch the error and try to compensate in some way -> adj.p
      q.vals <- try(qvalue::qvalue(p.vals)$qvalues, silent = TRUE)
      if (class(q.vals) == "try-error") {
        q.vals <- p.adjust(p.vals, "BH")
      }

      sig.genes <- rownames(exprs[-no.var.genes, , drop = F])[which(q.vals < 0.3)]
      # the genes are ordered by decreasing fold changes compared to other cell subsets
      temp.exprs <- exprs[-no.var.genes, , drop = F]
      fold.changes <- log2(rowMeans(temp.exprs[which(q.vals < 0.3), which(labs == 0), drop = F])) - log2(rowMeans(temp.exprs[which(q.vals < 0.3), which(labs == 1), drop = F]))
    } else {
      t.test.result <- multtest::mt.teststat(
        X = exprs,
        classlabel = labs,
        test = "t"
      )
      p.vals <- 2 * pt(abs(t.test.result), length(labs) - 2, lower.tail = FALSE)

      # catch same error as above
      q.vals <- try(qvalue::qvalue(p.vals)$qvalues, silent = TRUE)
      if (class(q.vals) == "try-error") {
        q.vals <- p.adjust(p.vals, "BH")
      }

      # calculate fold changes
      sig.genes <- rownames(exprs)[which(q.vals < 0.3)]
      genes.to.remove <- which(q.vals >= 0.3)
      if (length(genes.to.remove) > 0) {
        fold.changes <- log2(rowMeans(exprs[-genes.to.remove, which(labs == 0), drop = F])) - log2(rowMeans(exprs[-genes.to.remove, which(labs == 1), drop = F]))
      } else {
        fold.changes <- log2(rowMeans(exprs[, which(labs == 0), drop = F])) - log2(rowMeans(exprs[, which(labs == 1), drop = F]))
      }
    }
    # add genes for each type ordered by decreasing fold change
    deg.per.type[[t]] <- sig.genes[order(abs(fold.changes), decreasing = TRUE)]
  }

  # reduce to one reference profile per cell type
  ref.profiles <- matrix(
    NA,
    nrow = nrow(exprs),
    ncol = length(unique(pheno[, cell.type.column]))
  )
  colnames(ref.profiles) <- unique(pheno[, cell.type.column])
  rownames(ref.profiles) <- rownames(exprs)
  for (t in colnames(ref.profiles)) {
    ref.profiles[, t] <- rowMeans(
      exprs[, which(pheno[, cell.type.column] == t), drop = F]
    )
  }
  # again following Newman et al.:
  # take top g genes for every cell type, create signature matrices
  # choose the gene set that minimizes condition number
  limit <- min(max(sapply(deg.per.type, length)), max.genes)

  if (optimize) {
    cond.nums <- c()
    for (g in 1:limit) {
      all.genes <- c()
      for (sub.genes in deg.per.type) {
        all.genes <- c(all.genes, sub.genes[1:min(length(sub.genes), g)])
      }
      all.genes <- unique(all.genes)

      # several genes have 0 variance, therefore mt.multtest returns NA
      if (any(is.na(all.genes))) {
        all.genes <- all.genes[-which(is.na(all.genes))]
      }

      # estimate condition number of reduced matrix
      cond.nums <- c(cond.nums, kappa(ref.profiles[all.genes, ,drop=F], exact = FALSE))
    }
    optimal.g <- (1:limit)[which.min(cond.nums)]
  } else {
    optimal.g <- limit
  }

  # create gene list with the optimal g
  optimal.genes <- sapply(deg.per.type, function(x) {
    x[1:optimal.g]
  })
  optimal.genes <- unique(as.vector(optimal.genes))

  # deal with eventual NAs again
  if (any(is.na(optimal.genes))) {
    optimal.genes <- optimal.genes[-which(is.na(optimal.genes))]
  }

  # once again depending on whether split.data is true and possible
  if (split.data && !any(type.counts < 3)) {
    ref.mat <- sig.matrix[optimal.genes, ]
  } else {
    ref.mat <- ref.profiles[optimal.genes, ]
  }

  # remove any duplicate genes that might be in the matrix
  if(any(duplicated(rownames(ref.mat)))){
    warning("Found duplicates in reference matrix")
    ref.mat <- ref.mat[-which(duplicated(rownames(ref.mat))),]
  }

  return(sig.matrix = ref.mat)
}
