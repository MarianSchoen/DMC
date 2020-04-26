#' deconvolute given bulks with BSEQ-sc using single cell data
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
#' @example 
#' TODO: where do the 'training.exprs', 'training.pheno' and 
#' 'bulks.exprs' come from? 
#' 
#' run_bseqsc(
#'    training.exprs
#'    , training.pheno
#'    , bulks.exprs
#'  )
run_bseqsc <- function(
  exprs,
  pheno,
  bulks,
  exclude.from.signature = NULL,
  max.genes = NULL,
  optimize = TRUE,
  split.data = TRUE,
  verbose = FALSE,
  model = NULL
  ) {
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
  # ExpressionSet does strange things to names...
  rownames(exprs) <- make.names(rownames(exprs))
  rownames(bulks) <- make.names(rownames(bulks))
  
  bulk.pheno <- data.frame(colnames(bulks))
  rownames(bulk.pheno) <- colnames(bulks)
  colnames(bulk.pheno) <- "sample"
  colnames(bulks) <- rownames(bulk.pheno)
  bulks <- ExpressionSet(
    assayData = bulks,
    phenoData = AnnotatedDataFrame(bulk.pheno)
  )

  # do not normalize to counts for BSEQ-sc
  # exprs <- scale_to_count(exprs)

  if (!is.null(exclude.from.signature)) {
    to.exclude <- which(pheno[, "cell_type"] %in% exclude.from.signature)
    if  (length(to.exclude) > 0)  {
      exprs <- exprs[, -to.exclude]
      pheno <- pheno[-to.exclude, ]
    }
  }

  valid.model <- T
  if(!is.null(model)){
    if(all(c("ref.profiles", "g") %in% names(model))){
      full.mat <- model$ref.profiles
      g <- model$g
      if(all(names(g) %in% rownames(full.mat)) && length(g) == nrow(full.mat)){
        g <- g[rownames(full.mat)]
        B <- apply(full.mat, 2, function(x){x*g})
        # remove genes with g == 0 from reference
        if(any(rowSums(B) == 0)){
          B <- B[-which(rowSums(B) == 0),]
        }
      }else{
        warning("reference profiles and g vector do not contain the same genes")
        valid.model <- F
      }
    }else{
      warning("passed model parameter does not contain entries 'ref.profiles' and 'g'")
      valid.model <- F
    }
  }else{
    valid.model <- F
  }
  if(!valid.model){
    # use all supplied genes
    n.genes <- nrow(exprs)
    if (verbose) print("creating marker gene lists...")
    if (split.data) {
      learning.cells <- c()
      for (t in unique(pheno[, "cell_type"])) {
        learning.cells <- c(
          learning.cells,
          sample(
            which(pheno[, "cell_type"] == t),
            floor(0.9 * length(which(pheno[, "cell_type"] == t))),
            replace = FALSE
          )
        )
      }
  
      deg.per.type <- try(marker_genes(exprs[, learning.cells, drop = F], pheno[learning.cells, , drop = F], NULL))
      exprs <- exprs[, -learning.cells, drop = F]
      pheno <- pheno[-learning.cells, , drop = F]
    } else {
      deg.per.type <- try(marker_genes(exprs, pheno, NULL))
    }
    if (class(deg.per.type) == "try-error") {
      warning("Error while trying to find marker genes")
      return(list(est.props = NULL, sig.matrix = NULL, ref.profiles = NULL, g = NULL))
    }
  
    if (length(deg.per.type) == 0) {
      warning("No genes passed the BSEQ-sc criteria")
      return(list(est.props = NULL, sig.matrix = NULL, ref.profiles = NULL, g = NULL))
    }
  
    
    # create single cell and bulk expression set
    sc.exprs <- ExpressionSet(
      assayData = exprs,
      phenoData = AnnotatedDataFrame(pheno)
    )
  
    B <- try({
      bseqsc_basis(sc.exprs, deg.per.type,
        clusters = "cell_type", samples = "patient",
        ct.scale = TRUE
      )
    })
    
    all.per.type <- deg.per.type
    for(ct in names(all.per.type)){
      all.per.type[[ct]] <- rownames(exprs)
    }
    full.mat <- try({bseqsc_basis(sc.exprs, all.per.type,
                             clusters = "cell_type", samples = "patient",
                             ct.scale = TRUE)}
    )
    if (class(B) == "try-error") {
      warning("BSEQ-sc signature matrix creation failed")
      # save information in case of error
      #save(sc.exprs, deg.per.type, file = paste("../bseqsc-error_", Sys.time(), sep = ""))
      return(list(est.props = NULL, sig.matrix = NULL, ref.profiles = NULL, g = NULL))
    }
    g <- rep(0, nrow(full.mat))
    names(g) <- rownames(full.mat)
    g[rownames(B)] <- 1
  }
  fit <- try(bseqsc_proportions(bulks, B, log = F, verbose = verbose))
  if (class(fit) == "try-error") {
    warning("BSEQ-sc estimation failed")
    return(list(est.props = NULL, sig.matrix = NULL, ref.profiles = NULL, g = NULL))
  }
  if(!all(unique(pheno$cell_type) %in% rownames(fit$coefficients))){
    est.props <- complete_estimates(fit$coefficients, unique(pheno$cell_type))
  }else{
    est.props <- fit$coefficients
  }
  return(list(est.props = est.props, sig.matrix = B, ref.profiles = full.mat, g = g))
}
