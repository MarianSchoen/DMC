#' load results from temp_dir, create all possible plots and save them into
#' temp_dir/report_plots
#' 
#' @param temp_dir string, directory containing benchmark results
#' @param genesets list of cahracter vectors specifying gene sets used
#' @param features character vector containing all available features in data
#' @return NULL, save plots to disk
#' @export

plot_all <- function(
    temp_dir,
    genesets,
    features
  ) {
  plot_dir <- paste0(temp_dir, "/report_plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  # check for results on real bulks
  real_dir <- paste0(temp_dir, "/results/real/")
  real_files <- list.files(
    real_dir,
    full.names = TRUE,
    pattern = "deconv.*.h5"
  )
  if (length(real_files) > 0) {
    real_df <- data.frame()
    result_lists <- list()
    for (f in real_files) {
      result_lists[[f]] <- read_result_list(f)
      real_df <- rbind(
        real_df,
        prepare_data(result_lists[[f]])
      )
    }
  }
  
  bootstrap_files <- list.files(
    real_dir,
    full.names = TRUE,
    pattern = "bootstrap.*.h5"
  )
  if (length(bootstrap_files) > 0) {
    bootstrap_data <- NULL
    for (f in bootstrap_files) {
      bootstrap_data <- rbind(
        bootstrap_data,
        h5_read_mat(f)
      )
    }
    bootstrap_data <- as.data.frame(bootstrap_data)
    bootstrap_data$score <- as.numeric(as.character(bootstrap_data$score))
  }
  
  # create plots
  if (exists("real_df")) {
    
    # table plot
    score.plot.list <- evaluation_plot(
      real_df,
      "deconvolution quality"
    )
    score.plot <- score.plot.list$plot
    celltypes.ordered <- score.plot.list$celltype.order
    algorithms.ordered <- score.plot.list$algorithm.order
    
    # runtime
    runtime.plot <- plot_runtime(
      real_df, 
      algorithm.order = algorithms.ordered
    )
    
    # scatterplots
    scatter.plots.list <- list()
    for (i in seq_len(length(result_lists))) {
	    if("bulk.props" %in% names(result_lists[[i]])){
      scatter.plots.list[[i]] <- create_scatterplots(
        result_lists[[i]],
        celltype.order = celltypes.ordered,
        algorithm.order = algorithms.ordered
      )
	    }else{
		    cat("Something went wrong. bulk.props?")
	    	print(str(result_lists))
	    }
    }
    
    # condition number plots
    cond.num.plots <- plot_cond_num(
      real_df,
      algorithm.order = algorithms.ordered
    )
  }
  
  # bootstrap plots
  if (exists("bootstrap_data")) {
    bootstrap.plot <- ggplot2::ggplot(
      bootstrap_data, aes(x = algorithm, y = score, col = algorithm)
    ) + 
    facet_grid(rows = vars(cell_type)) +
    geom_boxplot() +
    ylim(0,1)
  }
  
  
  # simulation plots
  bulk_dir <- paste0(temp_dir, "/results/simulation/bulks")
  if (dir.exists(bulk_dir)) {
    results.lists <- list()
    files <- list.files(
      bulk_dir,
      full.names = TRUE,
      pattern = "*.h5"
    )
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        results.lists[[i]] <- read_result_list(files[i])
      }
      bulks.df <- data.frame()
      for (i in seq_len(length(results.lists))) {
        bulks.df <- rbind(
          bulks.df,
          prepare_data(results.lists[[i]])
        )
      }
      
      score.plot.sim.list <- evaluation_plot(
        bulks.df,
        "deconvolution quality (simulated bulks)"
      )
      
      score.plot.sim <- score.plot.sim.list$plot
      algorithms.ordered.sim <- score.plot.sim.list$algorithm.order
      celltypes.ordered.sim <- score.plot.sim.list$celltype.order
      
      runtime.plot.sim <- plot_runtime(
        bulks.df,
        algorithm.order = algorithms.ordered.sim
      )
      
      scatter.plots.list.sim <- list()
      for (i in seq_len(length(results.lists))) {
        scatter.plots.list.sim[[i]] <- create_scatterplots(
          results.lists[[i]],
          celltype.order = celltypes.ordered.sim,
          algorithm.order = algorithms.ordered.sim
        )
      }
    }
  }else{
    algorithms.ordered.sim <- NULL
    celltypes.ordered.sim <- NULL
  }
  
  # sample benchmark
  sample_dir <- paste0(temp_dir, "/results/simulation/samples")
  if (dir.exists(sample_dir)) {
    results.lists <- list()
    files <- list.files(
      sample_dir,
      full.names = TRUE,
      pattern = "*.h5"
    )
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        results.lists[[i]] <- read_result_list(files[i])
      }
      samples.df <- data.frame()
      for (i in seq_len(length(results.lists))) {
        samples.df <- rbind(
          samples.df,
          prepare_data(results.lists[[i]])
        )
      }
      sample.plots <- create_boxplots(
        samples.df,
        celltype.order = celltypes.ordered.sim,
        algorithm.order = algorithms.ordered.sim
      )
    }
  }
  
  # gene benchmark
  gene_dir <- paste0(temp_dir, "/results/simulation/genes")
  if (dir.exists(gene_dir)) {
    results.lists <- list()
    files <- list.files(
      gene_dir,
      full.names = TRUE,
      pattern = "*.h5"
    )
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        results.lists[[i]] <- read_result_list(files[i])
      }
      genes.df <- data.frame()
      for (i in seq_len(length(results.lists))) {
        genes.df <- rbind(
          genes.df,
          prepare_data(results.lists[[i]])
        )
      }
      gene.plots <- create_lineplots(
        genes.df,
        genesets,
        features,
        algorithm.order = algorithms.ordered.sim,
        celltype.order = celltypes.ordered.sim
      )
    }
  }
  
  # subtype benchmark
  subtype_dir <- paste0(temp_dir, "/results/simulation/subtypes")
  if (dir.exists(subtype_dir)) {
    results.lists <- list()
    files <- list.files(
      subtype_dir,
      full.names = TRUE,
      pattern = "*.h5"
    )
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        results.lists[[i]] <- read_result_list(files[i])
      }
      subtypes.df <- data.frame()
      for (i in seq_len(length(results.lists))) {
        subtypes.df <- rbind(
          subtypes.df,
          prepare_data(results.lists[[i]])
        )
      }
      score.plot.subtype <- subtype_plot(subtypes.df)
    }
  }

  # score benchmark
  score_dir <- paste0(temp_dir, "/results/simulation/scores")
  if (dir.exists(score_dir)) {
    results_file <- paste0(score_dir, "/deconv_output.h5")
    results.list <- rhdf5::h5dump(paste0(results_file))
    score_plots <- create_fit_plots(results.list)
    
    n_score_algorithms <- length(
	unique(unlist(sapply(results.list, function(x){unique(x$algorithm)})))
    )
    n_score_celltypes <- length(unique(results.list[[1]]$celltype))
  }
  
  # save plots
  if (exists("real_df")) {
    score_width <- max(2 * length(unique(real_df$cell_type)), 10)
    score_height <- max(2 * length(unique(real_df$algorithm)), 10)
    
    pdf(
      paste0(plot_dir, "/score_plot_real.pdf"),
      width = score_width,
      height = score_height
    )
    plot(score.plot)
    dev.off()
    
    pdf(
      paste0(plot_dir, "/scatter_plots_real.pdf"),
      width = score_width * 1.5,
      height = 6
    )
    for (l in scatter.plots.list) {
      if (!is.null(l)) {
        for (p in l) {
          if (!is.null(p)) plot(p)
        }
      }
    }
    dev.off()
    
    pdf(
      paste0(plot_dir, "/runtime_plot_real.pdf"),
      width = 16,
      height = 9
    )
    plot(runtime.plot)
    dev.off()
    
    pdf(
      paste0(plot_dir, "/condition_number_plots_real.pdf"),
      width = 16,
      height = 9
    )
    plot(cond.num.plots$cond_num_plot)
    plot(cond.num.plots$cond_vs_score)
    plot(cond.num.plots$variation_plot)
    dev.off()
    
    
    bootstrap.plot.width <- length(unique(bootstrap_data$algorithm)) * 2
    bootstrap.plot.height <- length(unique(bootstrap_data$cell_type)) * 2
    pdf(
      paste0(plot_dir, "/bootstrap_plot_real.pdf"),
      width = bootstrap.plot.width,
      height = bootstrap.plot.height
    )
    plot(bootstrap.plot)
    dev.off()
  }
  
  # bulk simulation
  if (exists("bulks.df")) {
    score.width <- max(length(unique(bulks.df$cell_type)) * 2, 10)
    score.height <- max(length(unique(bulks.df$algorithm)) * 2, 10)
    pdf(
      paste0(plot_dir, "/score_plot_simulated.pdf"),
      width = score.width,
      height = score.height
    )
    plot(score.plot.sim)
    dev.off()
    
    pdf(
      paste0(plot_dir, "/scatter_plots_simulated.pdf"),
      width = score.width * 1.5,
      height = 6
    )
    for (l in scatter.plots.list.sim) {
      for (p in l) {
        plot(p)
      }
    }
    dev.off()
  }
  
  # sample simulation
  if (exists("sample.plots")) {
    pdf(
      paste0(plot_dir, "/sample_plots_simulated.pdf"),
      width = 16,
      height = 9
    )
    for (p in sample.plots$cell.type.plots) plot(p)
    dev.off()
  }
  
  # geneset simulation
  if (exists("gene.plots")) {
    pdf(
      paste0(plot_dir, "/gene_time_plot_simulated.pdf"),
      width = 16,
      height = 9
    )
    plot(gene.plots$runtime.plot)
    dev.off()
  
  
    pdf(
      paste0(plot_dir, "/gene_score_plots_simulated.pdf"),
      width = 16,
      height = 9
    )
    for (p in gene.plots$cell.type.plots) plot(p)
    dev.off()
  }
  
  # subtype simulation
  if (exists("score.plot.subtype")) {
    pdf(paste0(plot_dir, "/subtype_plot.pdf"), width = 16, height = 8)
    plot(score.plot.subtype$overall)
    dev.off()
  }

  # score plots / fits
  if (exists("score_plots")) {
    height <- max(2 * n_score_algorithms, 10)
    width <- max(2 * n_score_celltypes, 10)
    pdf(
      paste0(plot_dir, "/correlation_fits.pdf"), 
      width = width, 
      height = height
    )
    
    for (p in score_plots$fits) {
      if (!is.null(p)) plot(p)
    }
    dev.off()

    pdf(
      paste0(plot_dir, "/fit_offsets.pdf"),
      width = width, height = height
    )
    plot(score_plots$offset_plot)
    dev.off()

   
    pdf(
      paste0(plot_dir, "/fit_errors.pdf"),
      width = width, height = height
    )
    plot(score_plots$error_plot)
    dev.off()

   
    pdf(
      paste0(plot_dir, "/fit_scores.pdf"),
      width = width, height = height
    )
    plot(score_plots$combined_plot)
    dev.off()
  }
  return(NULL)
}
