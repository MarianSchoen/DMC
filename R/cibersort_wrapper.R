# written by Tim Mirus

#' deconvolute given bulks with CIBERSORT using single cell data
#'
#' @param exprs matrix containing single cell profiles as columns
#' @param pheno phenotype data corresponding to the expression matrix.
#' Has to contain single cell labels in a column named 'cell_type'
#' @param bulks matrix containing bulk expression profiles as columns
#' @param exclude.from.signature list of cell types not to be included in the signature matrix
#' @param max.genes maximum number of genes that will be included in the signature for each celltype
#' @param optimize logical, should the signature matrix be optimized by condition number? If FALSE, max.genes genes will be used
#' @param split.data logical, should the training data be split for signature matrix creation? If TRUE, 10% of the data will be used to build
#' the signature matrix and the rest will be used to estimate the optimal features
#' @return list with one entry: est.props, matrix containing for each bulk the estimated fractions of the cell types contained
#' @example run_cibersort(training.exprs, training.pheno, bulk.exprs)

# source the CIBERSORT function (not available as package)
run_cibersort <- function(exprs,
                          pheno,
                          bulks,
                          exclude.from.signature = NULL,
                          max.genes = 500,
                          optimize = TRUE,
                          split.data = FALSE) {
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
    if (!is.null(max.genes) && max.genes == 0) {
        max.genes <- NULL
    }

    # normalize cell profiles to fixed counts
    exprs <- scale_to_count(exprs)
    # create signature matrix
    ref.profiles <- create_sig_matrix(
	exprs,
        pheno,
        exclude.from.signature,
        max.genes = max.genes,
        optimize = optimize,
        split.data = split.data
    )

    df.sig <- data.frame(GeneSymbol = rownames(ref.profiles))
    df.sig <- cbind(df.sig, ref.profiles)

    # CIBERSORT expects input to be supplied as .txt files
    dir.create("CIBERSORT/")
    write.table(df.sig,
        file = "CIBERSORT/signature_matrix.txt", quote = FALSE, row.names = FALSE,
        sep = "\t"
    )


    df.mix <- data.frame(GeneSymbol = rownames(bulks))
    df.mix <- cbind(df.mix, bulks)

    write.table(df.mix,
        file = "CIBERSORT/mixture.txt", quote = FALSE, row.names = FALSE,
        sep = "\t"
    )

    # call CIBERSORT; quantile normalization is recommended by the authors
    # switch off permutation, as I am not interested in p-values
    result <- CIBERSORT(
        sig_matrix = "CIBERSORT/signature_matrix.txt",
        mixture_file = "CIBERSORT/mixture.txt",
        QN = TRUE, perm = 0
    )

    # drop the additional information in the last 3 columns
    est.props <- t(result[1:ncol(bulks), -((ncol(result) - 2):ncol(result)), drop = FALSE])

    # CIBERSORT automatically stores the results in a file,
    # but we do not need it
    file.remove("CIBERSORT-Results.txt")
    file.remove("CIBERSORT/signature_matrix.txt")
    file.remove("CIBERSORT/mixture.txt")
    unlink("CIBERSORT/")

    return(list(est.props = est.props, sig.matrix = ref.profiles))
}
