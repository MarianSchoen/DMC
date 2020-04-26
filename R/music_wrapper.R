#' deconvolute given bulks with MuSiC using single cell data
#'
#' @param exprs non negative numeric matrix containing single cell profiles
#'  as columns and features as rows
#' @param pheno data.frame, with 'nrow(pheno)' must equal 'ncol(exprs)'. 
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature vector of strings of cell types not to be
#' included in the signature matrix
#' @param max.genes numeric, maximum number of genes that will be included in 
#' the signature for each celltype
#' @param optimize boolean, should the signature matrix be optimized by
#' condition number? If FALSE, max.genes genes will be used
#' @param split.data boolean, should the training data be split for signature
#' matrix creation? If TRUE, 10% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal
#' features
#' @param verbose boolean
#' @param model list containing two entries:
#' 1) ref.profiles - matrix containing reference profiles for all cell types in its columns
#' 2) g - weight vector for genes. For algorithms that do not assign weights to features,
#' this will consist of ones and zeroes, depending on wether a feature is included in the model or not
#' @return list with four entries: 
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm (features x cell types); can be calculated from ref.profiles and g
#' 3) ref.profiles - complete reference matrix (features x cell type); contains all genes unweighted
#' 4) g - named weight vector g; specifies for all genes, whether they are used in the effective signature (0,1) and
#' optionally assigns a weight to each gene (e.g. for DTD)
#' @example run_music(training.exprs, training.pheno, bulks[,1:5])
#'
run_music <- function(exprs,
                      pheno,
                      bulks,
                      exclude.from.signature = NULL,
                      max.genes = NULL,
                      optimize = TRUE,
                      split.data = TRUE,
                      model = NULL) {
  # error checking
  if (nrow(pheno) != ncol(exprs)) {
      stop("Number of columns in exprs and rows in pheno do not match")
  }
  if (nrow(exprs) != nrow(bulks)) {
      features <- intersect(rownames(exprs), rownames(bulks))
      if (length(features) > 0) {
          exprs <- exprs[features, ]
          bulks <- bulks[features, ]
      }
  }
  if(!is.null(max.genes)){
    if (max.genes == 0) {
        max.genes <- NULL
    }
  }


  # MuSiC uses all supplied genes
  n.genes <- nrow(exprs)

  # do not normalize profiles to counts
  # exprs <- scale_to_count(exprs)
  # ExpressionSet does strange things to names...

  rownames(exprs) <- make.names(rownames(exprs))
  rownames(bulks) <- make.names(rownames(bulks))
  # create ExpressionSets from exprs, pheno and bulks
  sc.exprs <- ExpressionSet(
    assayData = exprs,
    phenoData = AnnotatedDataFrame(pheno)
  )
  bulk.pheno <- data.frame(sample = colnames(bulks))

  # for some reason this is necessary, otherwise equality check of sample names will fail...
  rownames(bulk.pheno) <- colnames(bulks)
  colnames(bulks) <- rownames(bulk.pheno)
  bulks <- ExpressionSet(
    assayData = bulks,
    phenoData = AnnotatedDataFrame(bulk.pheno)
  )
  if (!is.null(exclude.from.signature)) {
    if (length(which(unique(pheno[, "cell_type"]) %in% exclude.from.signature)) > 0) {
      include.in.x <- unique(pheno[, "cell_type"])[-which(unique(pheno[, "cell_type"]) %in% exclude.from.signature)]
    }else{
      include.in.x <- unique(pheno[, "cell_type"])
    }
  } else {
    include.in.x <- unique(pheno[, "cell_type"])
  }

  # deconvolution
  est.prop.music <- try(music_prop(
    bulk.eset = bulks, sc.eset = sc.exprs, clusters = "cell_type",
    samples = "patient", select.ct = include.in.x, verbose = FALSE
  ))
  if(class(est.prop.music) == "try-error"){
	  return(list(est.props = NULL, sig.matrix = NULL, ref.profiles = NULL, g = NULL))
  }

  est.props <- as.matrix(t(est.prop.music$Est.prop.weighted))
  if(!all(include.in.x %in% rownames(est.props))){
    est.props <- complete_estimates(est.props, include.in.x)
  }
  full.mat <- NULL
  g <- NULL
  
  return(list(est.props = est.props, sig.matrix = NULL, ref.profiles = full.mat, g = g))
}
