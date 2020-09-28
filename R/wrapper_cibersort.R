#' deconvolute given bulks with CIBERSORT using single cell data
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
#' matrix creation? If TRUE, 10\% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal
#' features. default: FALSE
#' @param verbose boolean
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? 
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default NULL
#' @return list with four entries: 
#' 1) est.props - matrix containing for each bulk the
#' estimated fractions of the cell types contained
#' 2) sig.matrix - effective signature matrix used by the algorithm (features x cell types)

# source the CIBERSORT function (not available as package)
run_cibersort <- function(
    exprs,
    pheno,
    bulks,
    exclude.from.signature = NULL,
    max.genes = 500,
    optimize = TRUE,
    split.data = FALSE,
    cell.type.column = "cell_type",
    patient.column = NULL, 
    scale.cpm 
    ) {
	suppressMessages(library(e1071, quietly =TRUE))
	suppressMessages(library(parallel, quietly = TRUE))
	suppressMessages(library(preprocessCore, quietly = TRUE))
    # error checking
    if (nrow(pheno) != ncol(exprs)) {
        stop("Number of columns in exprs and rows in pheno do not match")
    }
    features <- intersect(rownames(exprs), rownames(bulks))
    if (length(features) > 0) {
        exprs <- exprs[features, ]
        bulks <- bulks[features, ]
    }
    if (!is.null(max.genes) && max.genes == 0) {
        max.genes <- NULL
    }

    
    if(scale.cpm){
        # prepare phenotype data and cell types to use
        exprs <- scale_to_count(exprs)
    }
    
    # create signature matrix
    ref.profiles <- create_sig_matrix(
	exprs,
        pheno,
        exclude.from.signature,
        max.genes = max.genes,
        optimize = optimize,
        split.data = split.data,
        cell.type.column = cell.type.column
    )
    df.sig <- data.frame(GeneSymbol = rownames(ref.profiles))
    df.sig <- cbind(df.sig, ref.profiles)

    # CIBERSORT expects input to be supplied as .txt files
    dir.create("CIBERSORT/")
    write.table(df.sig,
        file = "CIBERSORT/signature_matrix.txt", quote = FALSE, row.names = FALSE,
        sep = "\t"
    )

    # create data frame containing bulks
    df.mix <- data.frame(GeneSymbol = rownames(bulks))
    df.mix <- cbind(df.mix, bulks)
    write.table(df.mix,
        file = "CIBERSORT/mixture.txt", quote = FALSE, row.names = FALSE,
        sep = "\t"
    )

    # call CIBERSORT; quantile normalization is recommended by the authors
    # switch off permutation, as we are not interested in p-values
    result <- CIBERSORT(
        sig_matrix = "CIBERSORT/signature_matrix.txt",
        mixture_file = "CIBERSORT/mixture.txt",
        QN = TRUE, perm = 0
    )

    # drop the additional information in the last 3 columns
    est.props <- t(result[1:ncol(bulks), -((ncol(result) - 2):ncol(result)), drop = FALSE])

    # complete the estimation matrix in case of cell type dropouts
    if(!all(colnames(ref.profiles) %in% rownames(est.props))){
        est.props <- complete_estimates(est.props, colnames(ref.profiles))
    }
    # CIBERSORT automatically stores the results in a file,
    # but we do not need it
    file.remove("CIBERSORT-Results.txt")
    file.remove("CIBERSORT/signature_matrix.txt")
    file.remove("CIBERSORT/mixture.txt")
    unlink("CIBERSORT", recursive = TRUE)

    return(list(est.props = est.props, sig.matrix = ref.profiles))
}
