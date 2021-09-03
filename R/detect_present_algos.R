#' read results written by previous benchmarks and return a list
#' of algorithms already present in the data
#'
#' @param target_dir character, directory from which results are to be read
#' @param name_pattern character, pattern by which relevant files are recognized
#' @return character vector containing algorithm names

present_algos <- function(target_dir, name_pattern){
  previous_results <- list()
  if (dir.exists(target_dir)) {
    files <- list.files(
      target_dir,
      full.names = TRUE,
      pattern = name_pattern
    )
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        f <- files[i]
        previous_results[[i]] <- read_result_list(f)
      }
    }
  }else{
    dir.create(
      target_dir,
      recursive = TRUE
    )
  }

  if (length(previous_results) == 0) {
    return(NULL)
  }else{
    present_algorithms <- c()
    for (r in previous_results) {
      present_algorithms <- c(
        present_algorithms,
        unique(as.character(prepare_data(r)$algorithm))
      )
    }
    return(present_algorithms)
  }
}
