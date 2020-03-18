# written by Tim Mirus

split_dataset <- function(exprs, pheno, method = "random", prop = 0.25, grouping = NULL) {
    if(ncol(exprs) != nrow(pheno)){
        stop("expression and pheno data do not match")
    }
    # method can be 'random' or 'predefined'
    if (method == "random") {
        if(prop <= 0 | prop >= 1) {
            stop("Test set size out of bounds (0,1)")
        }
        test.samples <- sample(
        1:ncol(exprs),
            size = ceiling(0.25 * ncol(exprs)),
            replace = FALSE
        )
        training.samples <- (1:nrow(pheno))[-test.samples]
    }else if (method == "predefined") {
        if(length(grouping) != ncol(exprs) | length(unique(grouping)) != 2){
            stop("Please specify a valid grouping vector containing 1 and 2")
        }
        training.samples <- which(grouping == 1)
        test.samples <- which(grouping == 2)
    }else{
        stop("Please specify method to be either 'random' or 'predefined'")
    }
    # splitting
    test.exprs <- exprs[, test.samples]
    test.pheno <- pheno[test.samples, ]
    training.exprs <- exprs[, training.samples]
    training.pheno <- pheno[training.samples, ]

    return(list(training = list(exprs = training.exprs, pheno = training.pheno), 
                test = list(exprs = test.exprs, pheno = test.pheno))
            )
}
