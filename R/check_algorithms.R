#' check a list of algorithms for compatitbility with \link{benchmark}
#'
#' @param algorithms list containing a list for each algorithm. Each sublist
#' contains \cr 1) name: character  \cr 2) algorithm: function
#' @param feature_list, character vector of feature names for simulated data; default NULL
#' @return
#' @export

check_algorithms <- function(algorithms, feature_list = NULL) {
    if (!length(algorithms) > 0 || any(sapply(algorithms, length) < 3)) {
      stop("Check algorithm list. It seems to be incompatible.")
    }
    cat("Checking Algorithms for compatibility...\n")

    nFeatures <- 25
    if (!is.null(feature_list)) {
	    nFeatures <- length(feature_list)
    }

    # generate random scRNA-seq like data
    random.data <- DTD::generate_random_data(
        n.types = 3,
        n.samples.per.type = 60,
        n.features = nFeatures,
        sample.type = "Cell",
        feature.type = "gene",
        seed = 1234
    )
    rownames(random.data) <- feature_names
    # for each scRNA-seq profile, get its's "cell type"
    pheno.data <- sapply(
      strsplit(colnames(random.data), ".", fixed = TRUE),
      FUN = function(x) {
          x[[2]]
      }
    )
    names(pheno.data) <- colnames(random.data)
    patients <- sample(
      x = c("patient1", "patient2", "patient3"),
      size = 180,
      replace = TRUE
    )
    # cast the cell type information into a data.frame
    pheno <- data.frame(
      sample.name = names(pheno.data),
      cell_type = pheno.data,
      patient = patients
    )
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
    for (a in algorithms) {
      cat(a$name, "...\t")
      res <- try({
        a$algorithm(
          random.data,
          pheno,
          bulks$mixtures,
          cell.type.column = "cell_type",
          patient.column = "patient"
        )
      })

      if (length(class(res)) == 1) {
        if (class(res) == "try-error" || !is.list(res) ||
          !all(c("est.props", "sig.matrix", "model") %in% names(res))) {
            warning(
              paste(
                  "Algorithm ", a$name,
                  " did not return expected values
                  (est.props, sig.matrix, model). Please check implementation"
              )
            )
        }
      } else {
        warning(
          "Algorithm ", a$name, "did not return expected values. (est.props, sig.matrix, model). Please check implementation"
        )
      } 
    }
    cat("\n")
}
