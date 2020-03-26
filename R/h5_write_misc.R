#' write misc input
#'
#' @param genesets list of string vectors
#' @param algorithm.names 
#' @param function.call call object
#' @param grouping vector
#' @param filename string, where should the data be stored?
#'
#' @return NULL, function saves into ‘filename‘
write_misc_input <- function(genesets, algorithm.names, function.call, grouping, filename) {
  library(rhdf5)
  h5createFile(filename)
  if(!is.null(genesets)){
    h5createGroup(filename, "genesets")
    for(g in names(genesets)){
      h5write(as.vector(genesets[[g]]), filename, paste("genesets", g, sep = "/"))
    }
  }
  h5write(as.vector(algorithm.names), filename, "algorithms")
  h5write(as.vector(grouping), filename, "grouping")
  h5createGroup(filename, "function_call")
  function.args <- as.character(function.call)
  function.argnames <- names(as.list(function.call))
  h5write(as.vector(function.args), filename, "function_call/args")
  h5write(as.vector(function.argnames), filename, "function_call/argnames")
}

