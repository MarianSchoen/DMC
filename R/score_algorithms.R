#' fit correlation models to deconvolution results obtained on simulated bulks
#'
#' @param counts non-negative numeric matrix of single-cell counts
#' (genes x cells)
#' @param pheno data frame containing phenotype information corresponding
#' to the cells in `counts`
#' @param bulk_counts non-negative numeric matrix of single-cell counts
#' (genes x cells); profiles to be used for bulk creation
#' @param bulk_pheno data frame containing phenotype information corresponding
#' to the cells in `bulk_counts`
#' @param bulk_cell_types character vector containing names of cell types to
#' be included in simulated bulks; default NULL
#' @param exclude_from_signature character vector containing names of cell
#' types to be excluded from deconvolution models; default NULL
#' @param column_names list with 3 entries: \cr
#' 1) cell.type.column: string, name of column in pheno data frames containing
#' cell types\cr
#' 2) patient.column: string, name of column in pheno data frames containing
#' patient names / sample of origin of the single-cell profiles\cr
#' 3) sample.column: string, name of the column containing names of the
#' profiles contained in the count matrices (must correspond to
#' colnames(counts / bulk_counts))
#' @param algorithm_list input list of algorithms; see `benchmark` for info
#' @param nrep integer > 0, number of repetitions (for averaging and
#' error estimation)
#' @param nsets integer > 0, number of datasets to be drawn for deconvolution
#' @param nbulks integer > 0, number of bulks in each of the `nsets` data sets
#' @param datasets list of simulated bulk datasets as returned by this function;
#' optional, default NULL
#' @return list containing two elements:\cr
#' 1) results - data frame containing fit results
#' 2) datasets - list containing bulk datasets for all variances
#' @export

score_algorithms <- function(counts, pheno,
                             bulk_counts, bulk_pheno,
                             bulk_cell_types = NULL, 
                             exclude_from_signature = NULL,
                             column_names = list(
                               cell.type.column = "celltype",
                               patient.column = "patient",
                               sample.column = "sample"
                             ),
                             algorithm_list, nrep = 5, nsets = 4, nbulks = 250,
                             datasets = NULL) {
  # make sure counts and pheno data match
  if (!ncol(counts) == nrow(pheno)) {
    stop("Data dimensions do not match.")
  }
  
  # check column names
  for (cn in column_names) {
    if (!cn %in% colnames(pheno)) {
      stop(
        paste0("Column name '", cn, "' not found in pheno data frame.")
      )
    }
    if (!cn %in% colnames(bulk_pheno)) {
      stop(
        paste0("Column name '", cn, "' not found in pheno data frame.")
      )
    }
  }

  if (any(rownames(bulk_counts) != rownames(counts))) {
    stop("counts and bulk counts do not agree.")
  }

  # check algorithm list
  if (length(algorithm_list) < 1) {
    stop("No algorithms provided.")
  }

  # restrict to bulk cell types
  if (!is.null(bulk_cell_types)) {
    common_ct <- intersect(
      bulk_cell_types,
      unique(pheno[[column_names$cell.type.column]])
    )
    if (length(common_ct) > 0) {
      samples <- which(
        pheno[[column_names$cell.type.column]] %in% bulk_cell_types
      )
      counts <- counts[, samples]
      pheno <- pheno[samples, ]
    }
    if (length(common_ct) > 0) {
      samples <- which(
        bulk_pheno[[column_names$cell.type.column]] %in% bulk_cell_types
      )
      bulk_counts <- bulk_counts[, samples]
      bulk_pheno <- bulk_pheno[samples, ]
    }
  }

  # create bulks for several datasets
  # and create covariance matrix for each dataset
  variances <- as.integer(seq(from = 1, to = 100, length.out = 20))
  results <- list()
  C_all <- c()

  cat("Creating datasets\t", as.character(Sys.time()), "\n")
  
  if (is.null(datasets)) {
    datasets <- list()
    for (v in variances) {
      cat(v, "\t")
      datasets[[as.character(v)]] <- create_bulks(
        exprs = bulk_counts,
        pheno = bulk_pheno,
        cell.type.column = column_names$cell.type.column,
        n.bulks = 1000,
        n.profiles.per.bulk = 200,
        sum.to.count = TRUE,
        frac = v
      )

      if (is.null(C_all)) {
        C_all <- datasets[[as.character(v)]]$props
      } else {
        C_all <- cbind(
          C_all,
          datasets[[as.character(v)]]$props[rownames(C_all), ]
        )
      }
    }
  } else {
    for (v in variances) {
      if (is.null(C_all)) {
        C_all <- datasets[[as.character(v)]]$props
      } else {
        C_all <- cbind(
          C_all,
          datasets[[as.character(v)]]$props[rownames(C_all), ]
        )
      }
    }
  }
  cat("\n")

  # create an overall covariance matrix
  results[["covariances"]] <- create_cov_mat(C_all)
  rm(C_all)

  # deconvolve each dataset with all algorithms
  cat("Start deconvolution\t", as.character(Sys.time()), "\n")
  for (v in names(datasets)) {
    cat(v, ":\t")
    results[[v]] <- list()
    for (a in names(algorithm_list)) {
      cat(a, "...\t")
      results[[v]][[a]] <- algorithm_list[[a]]$algorithm(
        exprs = counts,
        pheno = pheno,
        bulks = datasets[[v]]$bulks,
        cell.type.column = column_names$cell.type.column,
        exclude.from.signature = exclude_from_signature,
        patient.column = column_names$patient.column
      )$est.props
    }
    cat("\n")
  }
  deconv_results <- results

  cat(
    "Fit parameters to correlation curves...\t", as.character(Sys.time()), "\n"
  )

  # fit
  results <- data.frame()
  for (i in seq_len(nrep)) {
    fit_results <- fit_data(
      deconv_results,
      datasets,
      nbulks = nbulks,
      nsets = nsets
    )
    results <- rbind(
      results,
      cbind(fit_results, repetition = i)
    )
  }
  return(list(results = results, datasets = datasets))
}
