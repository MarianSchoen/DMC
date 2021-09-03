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
#' max.genes = NULL and all genes will be selected as maximum
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? default "cell_type"
#' @return numeric matrix, holding reference profiles in its column,
#' features in its rows

create_sig_matrix <- function(
  exprs,
  pheno,
  exclude.celltypes = NULL,
  max.genes = NULL,
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
    if (length(to.exclude) > 0) {
      exprs <- exprs[, -to.exclude, drop = F]
      pheno <- pheno[-to.exclude, , drop = F]
    }
  }

  # make sure that not more genes than available are selected
  if (is.null(max.genes)) {
	   max.genes <- nrow(exprs)
  }
  cat("Maximum number of genes per cell type: ", max.genes, "\n")

  # for each cell type test against all others for DEG
  # using two-sided t-test
  deg.per.type <- list()
  for (t in unique(pheno[, cell.type.column])) {
    labs <- ifelse(pheno[, cell.type.column] == t, 0, 1)
    t.test.result <- multtest::mt.teststat(
      X = exprs,
      classlabel = labs,
      test = "t"
    )
    p.vals <- 2 * pt(abs(t.test.result), length(labs) - 2, lower.tail = FALSE)
    names(p.vals) <- rownames(exprs)

    nna.genes <- names(p.vals)[!is.na(p.vals)]
    p.vals <- p.vals[nna.genes]

    # control FDR either via q-value or with benjamini-hochberg
    q.vals <- try({
      qvalue::qvalue(p.vals)$qvalues
    })
    if (class(q.vals) == "try-error") {
      warning("q-value calculation failed. Using BH-correction.")
    	q.vals <- p.adjust(p.vals, "BH")
      sig.entries <- which(q.vals < 0.1)
    }else{
      sig.entries <- which(q.vals < 0.3)
    }
    sig.genes <- nna.genes[sig.entries]

    # catch possible errors related to sig.genes
    if (any(is.na(sig.genes))) {
      sig.entries <- sig.entries[!is.na(sig.genes)]
      sig.genes <- sig.genes[!is.na(sig.genes)]
    }
    if (!length(sig.genes) > 0) break

    # calculate fold changes
    fold.changes <- log2(Matrix::rowMeans(
        exprs[sig.entries, which(labs == 0), drop = F]
      )) - log2(Matrix::rowMeans(
        exprs[sig.entries, which(labs == 1), drop = F]
      ))
    if (any(is.infinite(fold.changes))) {
      sig.genes <- sig.genes[-which(is.infinite(fold.changes))]
      fold.changes <- fold.changes[-which(is.infinite(fold.changes))]
    }
    if (any(is.nan(fold.changes))) {
      sig.genes <- sig.genes[-which(is.nan(fold.changes))]
      fold.changes <- fold.changes[-which(is.nan(fold.changes))]
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
  if (length(deg.per.type) > 0) {
  limit <- min(
    max(sapply(deg.per.type, length)), max.genes
  )
  }else{
    warning("No significant genes found. Returning NULL.")
    return(NULL)
  }

  cond.nums <- rep(Inf, times = limit)
  # condition number is unstable for very few genes, therefore use at least 10
  for (g in 10:limit) {
    all.genes <- unique(unlist(
      sapply(deg.per.type, function(sub.genes, lim = g) {
        sub.genes[1:min(length(sub.genes), lim)]
      })
    ))
    if (any(is.na(all.genes))) {
      all.genes <- all.genes[-which(is.na(all.genes))]
    }

    # estimate condition number of reduced matrix
    cond.nums[g] <- kappa(ref.profiles[all.genes,, drop = FALSE], exact = FALSE)
  }
  optimal.g <- which.min(cond.nums)

  # alternative:
  #optimal.g <- limit
  #cat(
  #  "Chose ", optimal.g,
  #  "genes per cell type resulting in condition number of",
  #  cond.nums[optimal.g], "\n"
  #)

  # create gene list with the optimal g
  optimal.genes <- unique(unlist(
    sapply(deg.per.type, function(sub.genes, opt.g = optimal.g) {
      sub.genes[1:min(opt.g, length(sub.genes))]
    })
  ))

  ref.profiles <- ref.profiles[optimal.genes, ]

  # remove any duplicate genes that might be in the matrix
  if (any(duplicated(rownames(ref.profiles)))) {
    warning("Found duplicates in reference matrix")
    ref.profiles <- ref.profiles[-which(duplicated(rownames(ref.profiles))), ]
  }
  return(ref.profiles)
}
