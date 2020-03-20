check_algorithms <- function(algorithms){
    randomRNA <- DTD::generate_random_data(
        n.types = 3,
        n.samples.per.type = 10,
        n.features = 50,
        sample.type = "Cell",
        feature.type = "gene",
        seed = 1234
    )

    pheno.data <- sapply(strsplit(colnames(randomRNA), ".", fixed = TRUE),
    FUN = function(x) {
        x[[2]]
    }
    )
    names(pheno.data) <- colnames(randomRNA)
    pheno <- data.frame(sample.name = names(pheno.data), cell_type = pheno.data)
    rownames(pheno) <- pheno$sample.name

    # create bulks from random data
    bulks <- DTD::mix_samples(
        expr.data = randomRNA,
        pheno = pheno.data,
        included.in.X = unique(pheno.data),
        n.samples = 50,
        n.per.mixture = 5
    )

    for(a in algorithms){
        res <- a$algorithm(
            randomRNA,
            pheno,
            bulks$mixtures
        )
        if(!is.list(res)){
            if(!all(c("est.props", "sig.matrix") %in% names(res))){
                stop("Algorithm did not return expected values. Please check implementation")
            }
        }
    }
}