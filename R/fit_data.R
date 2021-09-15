#' fit the correlation model of DAB to deconvolution results
#' 
#' @param results List containing one sublist for bulk dataset 
#' (each simulated with a different variance in cellular composition). 
#' Each sublist contains deconvolution results (matrix, cell type x bulk) for
#' all algorithms.
#' @param dataset The bulk datasets that were deconvoluted. List containing 
#' bulk counts and proportions as created by `create_bulks` for different 
#' variances
#' @param nsets integer > 0, number of datasets to be drawn for deconvolution
#' @param nbulks integer > 0, number of bulks in each of the `nsets` data sets
#' @return data frame containing results 
#' (correlations, model parameters, predictions ...)

fit_data <- function(results, datasets, nbulks = NULL, nsets = NULL) {
  algo_params <- list()
  params <- c(0.05, 0.05)
  names(params) <- c("E", "offset")
  for (algo in names(results[["1"]])) {
    algo_params[[algo]] <- list()
    cts <- rownames(results[["1"]][[algo]])
    for (ct in cts) {
      algo_params[[algo]][[ct]] <- params
    }
  }


  # prepare data for fitting
  df <- c()
  for (i in names(results)) {
    if (i == "covariances") next
    for (j in 1:nsets) {
      # create subset of bulks
      bulks <- sample(
        1:ncol(datasets[[i]]$props),
        size = nbulks, replace = TRUE
      )

      real_props <- datasets[[i]]$props[, bulks]
      # covariances
      cov_mat <- create_cov_mat(real_props)

      for (algo in names(results[[i]])) {
        # estimate for each algorithm
        est_props <- results[[i]][[algo]][, bulks]
        cts <- intersect(
          rownames(real_props),
          rownames(results[[i]][[algo]])
        )

        for (ct in cts) {
          ct_var <- cov_mat[ct, ct]
          score <- cor(real_props[ct, ], est_props[ct, ])

          info <- c(ct, algo, i, score, ct_var, j)
          df <- rbind(df, info)
        }
      }
    }
  }
  colnames(df) <- c(
    "celltype", "algorithm", "index", "correlation",
    "variance", "set"
  )
  df <- as.data.frame(df)

  # adjust variable types
  df$index <- as.numeric(df$index)
  df$correlation <- as.numeric(df$correlation)
  df$variance <- as.numeric(df$variance)

  # estimate standard deviation
  estimate_sd <- TRUE
  if (estimate_sd) {
    sds <- c()
    # estimate sd for each data point
    for (i in 1:nrow(df)) {
      v <- df[i, "variance"]
      ct <- df[i, "celltype"]
      al <- df[i, "algorithm"]
      sub_df <- df[intersect(which(df$algorithm == al), which(df$celltype == ct)), ]
      sorted_vars <- sort(sub_df$variance)
      # select data points with similar underlying variances
      var_pos <- which(sorted_vars == v)
      variances <- sorted_vars[max(0, var_pos - 3):min(var_pos + 3, length(sorted_vars))]
      # calculate sd from correlations in the bin
      sds <- c(sds, sd(sub_df[which(sub_df$variance %in% variances), "correlation"]))
    }
    df$error <- sds
  }

  # minimization
  algo_results <- list()

  for (algo in unique(df$algorithm)) {
    algo_results[[algo]] <- list()
    cat(algo, "...\n\t")
    algo_df <- df[which(df$algorithm == algo), ]



    temp_df <- algo_df
    offsets_algo <- c()
    E_algo <- c()
    for (ct in unique(temp_df$celltype)) {
      cat(ct, "...\t")
      opt_result <- optim(
        algo_params[[algo]][[ct]],
        fn = cost_fun,
        method = "Nelder-Mead",
        control = list(maxit = 1000),
        temp_df = temp_df,
        ct = ct
      )


      offsets_algo <- c(offsets_algo, abs(opt_result$par["offset"]))
      E_algo <- c(E_algo, abs(opt_result$par["E"]))
    }
    names(E_algo) <- unique(temp_df$celltype)
    names(offsets_algo) <- unique(temp_df$celltype)
    algo_results[[algo]] <- list(offsets = offsets_algo, errors = E_algo)
    cat("\n")
  }

  # predict values based on model
  fit_df <- data.frame()
  for (algo in unique(df$algorithm)) {
    temp_df <- df[which(df$algorithm == algo), ]
    for (ct in unique(temp_df$celltype)) {
      sub_df <- temp_df[which(temp_df$celltype == ct), ]
      predictions <- predict_value(
        sub_df, algo_results[[algo]]$offsets[ct],
        algo_results[[algo]]$errors[ct]
      )
      fit_df <- rbind(
        fit_df,
        data.frame(
          correlation = sub_df$correlation,
          prediction = predictions,
          index = sub_df$index,
          variance = sub_df$variance,
          algorithm = rep(algo, length(predictions)),
          celltype = rep(ct, length(predictions)),
          offset = algo_results[[algo]]$offsets[[ct]],
          error_sd = algo_results[[algo]]$errors[[ct]]
        )
      )
    }
  }

  # return raw data
  return(fit_df)
}
