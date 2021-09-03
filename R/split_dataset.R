#' split dataset in training and test set either by grouping vector or randomly
#'
#' @param exprs non-negative numeric matrix with features as rows, and
#' scRNA-Seq profiles as columns. 'ncol(sc.counts)' must equal 'nrow(sc.pheno)'
#' @param pheno data frame with scRNA-Seq profiles as rows, and pheno entries
#'  in columns. 'nrow(sc.pheno)' must equal 'ncol(sc.counts)'
#' @param method string specifying method of splitting;
#' either 'random' (default) or 'predefined'
#' @param prop numerical, fraction of profiles to be used as test set
#' @param grouping factor with 2 levels (training == 1, test == 2),
#' and 'length(grouping)' must be 'ncol(sc.counts)'.
#'  Assigns each scRNA-Seq profile to either test or train cohort.
#' @return list containing two lists:
#' 1) training
#' 2) test
#' both lists contain two entries:
#' 1) exprs - non-negative numeric matrix with features as rows, and
#' scRNA-Seq profiles as columns. 'ncol(sc.counts)' must equal 'nrow(sc.pheno)'
#' 2) pheno - data frame with scRNA-Seq profiles as rows, and pheno entries
#'  in columns. 'nrow(sc.pheno)' must equal 'ncol(sc.counts)'

split_dataset <- function(
  exprs,
  pheno,
  method = "random",
  prop = 0.25,
  grouping = NULL
) {
  # parameter check
  if (ncol(exprs) != nrow(pheno)) {
      stop("expression and pheno data do not match")
  }
  # method can be 'random' or 'predefined'
  if (method == "random") {
      # randomly assign given fraction of the samples to the test set
      if (prop <= 0 | prop >= 1) {
          stop("Test set size out of bounds (0,1)")
      }
      test.samples <- sample(
        seq_len(ncol(exprs)),
        size = ceiling(0.25 * ncol(exprs)),
        replace = FALSE
      )
      training.samples <- (seq_len(nrow(pheno)))[-test.samples]
  }else if (method == "predefined") {
      # divide dataset according to the grouping vector
      if (length(grouping) != ncol(exprs) | length(unique(grouping)) != 2 |
         !all(c(1,2) %in% grouping)) {
          stop("Please specify a valid grouping vector
                containing 1 (training) and 2 (test)")
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

  return(list(
    training = list(exprs = training.exprs, pheno = training.pheno),
    test = list(exprs = test.exprs, pheno = test.pheno)
  ))
}
