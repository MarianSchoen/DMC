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
#' @example run_cibersort(training.exprs, training.pheno, bulk.exprs)

# source the CIBERSORT function (not available as package)
run_cibersort <- function(exprs,
                          pheno,
                          bulks,
                          exclude.from.signature = NULL,
                          max.genes = 500,
                          optimize = TRUE,
                          split.data = FALSE,
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
    if (!is.null(max.genes) && max.genes == 0) {
        max.genes <- NULL
    }

    # normalize cell profiles to fixed counts
    exprs <- scale_to_count(exprs)
    # create signature matrix
    valid.model <- T
    if(!is.null(model)){
        if(all(c("ref.profiles", "g") %in% names(model))){
            full.mat <- model$ref.profiles
            g <- model$g
            if(all(names(g) %in% rownames(full.mat)) && length(g) == nrow(full.mat)){
                g <- g[rownames(full.mat)]
                if(any(duplicated(rownames(full.mat)))){
                    g <- g[-which(duplicated(rownames(full.mat)))]
                    full.mat <- full.mat[-which(duplicated(rownames(full.mat))),]
                }
                ref.profiles <- apply(full.mat, 2, function(x){x*g})
                if(any(rowSums(ref.profiles) == 0)){
                    ref.profiles <- ref.profiles[-which(rowSums(ref.profiles) == 0),]
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
        model <- create_sig_matrix(
    	exprs,
            pheno,
            exclude.from.signature,
            max.genes = max.genes,
            optimize = optimize,
            split.data = split.data
        )
        ref.profiles <- model$sig.matrix
        full.mat <- model$full.matrix
        g <- model$g
    }

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
    if(!all(colnames(ref.profiles) %in% rownames(est.props))){
        est.props <- complete_estimates(est.props, colnames(ref.profiles))
    }
    # CIBERSORT automatically stores the results in a file,
    # but we do not need it
    file.remove("CIBERSORT-Results.txt")
    file.remove("CIBERSORT/signature_matrix.txt")
    file.remove("CIBERSORT/mixture.txt")
    unlink("CIBERSORT", recursive = TRUE)

    return(list(est.props = est.props, sig.matrix = ref.profiles, ref.profiles = full.mat, g = g))
}
