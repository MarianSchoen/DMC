#' fine and coarse subtype benchmark function
#'
#' # basically, there are 3 different things to compare.
#' 1. correlation per algorithm on the deep C
#' 2. correlation per algorithm on the accumulated cursory C
#' 3. special scenario where subtype 1 is in X, and subtype 2 is in Y
#'
#' 1. and 2. can be done in one go, 3. needs more prepration.
#'
#' @param sc.counts count matrix, features as rows,
#' scRNA-Seq profiles as columns
#' @param sc.pheno data.frame. scRNA-Seq profiles as rows. Must have
#' 'cell.type.column' and 'sample.name.column'
#' @param cell.type.column string, column of 'sc.pheno' holding the
#' input cell type labels. Within these entries, the clustering is done.
#' @param algorithm.list List containing a list for each algorithm.
#' Each sublist contains 1) name, 2) function and 3) model
#' @param subtype.pattern character,
#' string by which subtype column is recognized; default "subtype"
#' @param sample.name.column string, column of the 'colnames(sc.counts)'
#' @param n.clusters integer vector of clustering depths (number of subclusters
#' created for each cell type), default c(1, 2, 4, 8).
#' This means that in the finest clustering, each celltype will be split in
#' 8 subtypes, in the next step each will be split in 4 subtypes, ...
#' @param verbose logical, should information about the process be printed?
#' @param patient.column string, which column of 'pheno'
#' holds the patient information; optional, default "patient"
#' @param n.bulks numeric > 0, number of bulks to simulate. default 500
#'
#' @return list containing deconvolution results for
#' different cell type granularities

fine_coarse_subtype_benchmark <- function(
  sc.counts,
  sc.pheno,
  algorithm.list,
  subtype.pattern = "subtype",
  cell.type.column = "cell_type",
  sample.name.column = "sample.name",
  n.clusters = c(1, 2, 4, 8),
  verbose = FALSE,
  patient.column = "patient",
  n.bulks = 500
) {
  if (ncol(sc.counts) != nrow(sc.pheno)) {
    stop(
      "In DAB::fine_coarse_subtype_benchmark: 'sc.counts' does not fit
      'sc.pheno'. 'ncol(sc.counts)' must equal 'nrow(sc.pheno)'"
    )
  }

  if (! cell.type.column %in% colnames(sc.pheno)) {
    stop(
      "In DAB::fine_coarse_subtype_benchmark: 'cell.type.column' is not in
      'colnames(sc.pheno)'"
    )
  }

  if (! sample.name.column %in% colnames(sc.pheno)) {
    stop(
      "In DAB::fine_coarse_subtype_benchmark: 'sample.name.column' is not in
      'colnames(sc.pheno)'"
    )
  }

  # find a colname for the subtype information:
  while (any(grepl(pattern = subtype.pattern, x = colnames(sc.pheno)))) {
    subtype.pattern <- paste0(subtype.pattern, sample.int(100, 1))
  }

  # subcluster the data, return the sc.pheno frame with additional columns
  sc.pheno <- add_subtype_pheno(
    sc.counts = sc.counts,
    sc.pheno = sc.pheno,
    cell.type.column = cell.type.column,
    sample.name.column = sample.name.column,
    new.subtype.column = subtype.pattern,
    hclust.obj = NA,
    n.clusters = n.clusters,
    verbose = verbose
  )

  # detect the "new" subtype columns:
  subtype.columns <- colnames(sc.pheno)[
    grepl(pattern = subtype.pattern, x = colnames(sc.pheno))
  ]

  # initialise the output list.
  # (this holds as many entries as there are 'subtype.columns')
  # each entry holds the following entries:
  #  1. a list entry for each 'subtype.column'
  #  2. There are four entries:
  #   - 'c.true': matrix, true cell compositions, for all subtypes
  #   - 'c.true.coarsly': matrix, true cell compositions, for the
  #                       'cell.type.column'cell types
  #   - 'c.estimated.list': a list, with a matrix entry for each algorithm.
  #          There, the estimated cell compositions, for all subtypes are stored
  #   - 'c.estimated.coarsly.list': a list, with a matrix entry per algorithm.
  #          There, the estimated cell compositions, for the
  #          'cell.type.column'cell types are stored.
  column.list <- vector(
    mode = "list",
    length = length(subtype.columns)
  )
  names(column.list) <- subtype.columns

  # go through each 'subtype.column'
  if (verbose) cat("Deconvoluting..\n")
  for (column in subtype.columns) {
    # simulate bulks, and deconvolute them
    if (verbose) cat(column, "\n")
    some.estimates <- deconvolute(
      training.expr = sc.counts
      , training.pheno = sc.pheno
      , test.expr = sc.counts
      , test.pheno = sc.pheno
      , cell.type.column = column
      , algorithms = algorithm.list
      , verbose = TRUE
      , n.repeats = 1
      , patient.column = patient.column
      , n.bulks = n.bulks
    )

    # extract the true c (with all subtype quantities)...
    c.true <- some.estimates$bulk.props
    # this C is quite deep. I need the coarse C, therefore:
    c.true.coarsly <- matrix(
      data = NA
      , nrow = length(unique(sc.pheno[, cell.type.column]))
      , ncol = ncol(c.true)
    )
    colnames(c.true.coarsly) <- colnames(c.true)
    rownames(c.true.coarsly) <- unique(sc.pheno[, cell.type.column])
    for (major.cell.type in rownames(c.true.coarsly)) {
      associated.subtypes <- rownames(c.true)[
        which(grepl(pattern = major.cell.type, x = rownames(c.true)))
      ]
      if (length(associated.subtypes) > 0) {
        c.true.coarsly[major.cell.type, ] <- colSums(
          x = c.true[associated.subtypes, , drop = FALSE]
        )
      }else{
        c.true.coarsly[major.cell.type, ] <- 0
      }
    }

    # store the true C matrices:
    column.list[[column]][["c.true"]] <- c.true
    column.list[[column]][["c.true.coarsly"]] <- c.true.coarsly

    # TODO: manage the "n.repeats = 1" problem
    algorithms <- names(some.estimates$results.list$`1`)

    # for the estimated C, initialise two list entries:
    column.list[[column]][["c.estimated.list"]] <- vector(
      mode = "list"
      , length = length(algorithms)
    )
    names(column.list[[column]][["c.estimated.list"]]) <- algorithms


    column.list[[column]][["c.estimated.coarsly.list"]] <- vector(
      mode = "list"
      , length = length(algorithms)
    )
    names(column.list[[column]][["c.estimated.coarsly.list"]]) <- algorithms

    # go through all provided algorithms
    for (algorithm in algorithms) {
      # for the current algorithm, extract the estimated 'C's
      c.estimated <- some.estimates$results.list$`1`[[algorithm]]$est.props
      # store it:
      column.list[[column]][["c.estimated.list"]][[algorithm]] <- c.estimated

      # add up those subtypes, that origin from the same cell type:
      c.estimated.coarsly <- matrix(
        data = NA
        , nrow = nrow(c.true.coarsly)
        , ncol = ncol(c.true.coarsly)
      )
      dimnames(c.estimated.coarsly) <- dimnames(c.true.coarsly)
      for (major.cell.type in rownames(c.true.coarsly)) {
        associated.subtypes <- rownames(c.estimated)[
          which(grepl(pattern = major.cell.type, x = rownames(c.estimated)))
        ]
        if (length(associated.subtypes) > 0) {
          c.estimated.coarsly[major.cell.type, ] <- colSums(
            x = c.estimated[associated.subtypes, , drop = FALSE]
          )
        }else{
          c.estimated.coarsly[major.cell.type, ] <- 0
        }
      }
      # and store again
      column.list[[column]][["c.estimated.coarsly.list"]][[
        algorithm
      ]] <- c.estimated.coarsly
    }
  }
  return(column.list)
}
