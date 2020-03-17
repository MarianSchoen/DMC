# written by Tim Mirus

#' simulate bulk data by summing over single-cell data
#'
#' @param exprs expression matrix of the single cell data,
#' one single cell per column
#' @param pheno DataFrame, phenotype data containing
#' cell type labels for the expression matrix,
#' label column must be named 'cell_type',
#' ordering of cells must be the same as in exprs
#' @param n.bulks integer, the number of bulks to be created, defaults to 500
#' @param include.in.bulks list of cell types to be used for bulk simulation;
#' if not supplied, all will be used
#' @param fraction.per.bulk fraction of samples to be randomly
#' drawn for each bulk; default 0.1
#' @param sum.to.count logical, should all bulks be normalized
#' to a fixed total count number? default TRUE
#' @return list with two entries:
#' 1) matrix containing bulk expression profiles (features x bulks)
#' 2) matrix containing quantities (cell type x bulks)
#' @example create_bulks(training.exprs, training.pheno, n.bulks = 1000)

create_bulks <- function(exprs, pheno, n.bulks = 500, include.in.bulks = NULL, fraction.per.bulk = 0.1, sum.to.count = TRUE) {
    # error checking
    if (nrow(pheno) != ncol(exprs)) {
        stop("Number of columns in exprs and rows in pheno do not match")
    }

    # keep only specified cell types
    if (is.null(include.in.bulks)) {
        include.in.bulks <- unique(pheno[, "cell_type"])
    }

    if (length(which(pheno[, "cell_type"] %in% include.in.bulks)) > 0) {
        pheno <- pheno[which(pheno[, "cell_type"] %in% include.in.bulks), , drop = F]
    }
    exprs <- exprs[, rownames(pheno), drop = F]

    exprs <- scale_to_count(exprs)

    genes <- rownames(exprs)

    bulk.exprs <- matrix(
        0,
        nrow = nrow(exprs),
        ncol = n.bulks
    )

    rownames(bulk.exprs) <- genes
    colnames(bulk.exprs) <- as.character(1:n.bulks)

    # create a matrix to contain true proportions for each simulated bulk
    props <- matrix(0, nrow = length(unique(pheno[, 4])), ncol = n.bulks)
    rownames(props) <- unique(pheno[, 4])
    colnames(props) <- colnames(bulk.exprs)

    # create random bulks
    for (i in 1:ncol(bulk.exprs)) {
        bulk.samples <- sample(
            1:ncol(exprs),
            ceiling(fraction.per.bulk * ncol(exprs)),
            replace = TRUE
        )
        bulk.exprs[, i] <- rowSums(exprs[, bulk.samples, drop = F])
        # store the quantities
        for (t in rownames(props)) {
            props[t, i] <- sum(pheno[bulk.samples, "cell_type"] == t) / length(bulk.samples)
        }
    }
    if(any(duplicated(rownames(bulk.exprs)))){
        print("Found duplicate features in simulated bulks. Removing...")
        bulk.exprs <- bulk.exprs[-which(duplicated(bulk.exprs)), ]
    }

    if (sum.to.count) {
        bulk.exprs <- apply(bulk.exprs, 2, function(x) {
            x / sum(x) * length(x)
        })
    }

    return(list(bulks = bulk.exprs, props = props))
}
