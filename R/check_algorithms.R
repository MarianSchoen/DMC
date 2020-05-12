#' check_algorithms
#'
#' @param algorithms list containing a list for each algorithm. Each sublist
#' contains 1) name  and 2) function 
#' @return
#' @export
#' 
check_algorithms <- function(algorithms){
    if (!length(algorithms) > 0 || any(sapply(algorithms, length) != 2)) {
    stop("Check algorithm list")
  }
    # generate random scRNA-seq like data 
    random.data <- DTD::generate_random_data(
        n.types = 3,
        n.samples.per.type = 20,
        n.features = 25,
        sample.type = "Cell",
        feature.type = "gene",
        seed = 1234
    )
    # for each scRNA-seq profile, get its's "cell type"
    pheno.data <- sapply(strsplit(colnames(random.data), ".", fixed = TRUE),
    FUN = function(x) {
        x[[2]]
    }
    )
    names(pheno.data) <- colnames(random.data)
    # cast the cell type information into a data.frame
    pheno <- data.frame(sample.name = names(pheno.data), cell_type = pheno.data)
    rownames(pheno) <- pheno$sample.name

    # create bulks from random data
    bulks <- DTD::mix_samples(
        expr.data = random.data,
        pheno = pheno.data,
        included.in.X = unique(pheno.data),
        n.samples = 20,
        n.per.mixture = 4
    )

    # loop over all algorithms, run them on the random data, 
    # and check wheter the output is valid
    for(a in algorithms){
        res <- a$algorithm(
            random.data,
            pheno,
            bulks$mixtures
        )
    	
        if(!is.list(res) || !all(c("est.props", "sig.matrix") %in% names(res))){
            stop(
                paste(
                    "Algorithm ", a$name, " did not return expected values (est.props, sig.matrix). 
                    Please check implementation"
                    )
            )
        }
    }
}
