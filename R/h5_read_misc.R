#' read_misc_input
#'
#' @param filename string, which file should be read in?
#'
#' @return list with: 
#'    - 'genesets': list of string vectors
#'    - 'algorithms': vector containing algorithm names
#'    - 'grouping': numeric vector grouping samples in test and training set
#'    - 'function.call': call object

read_misc_input <- function(filename){
  # parameter check
  if(!file.exists(filename)) {
    stop(paste("Could not find file ", filename, sep = ""))
  }

  # read the content
  content <- rhdf5::h5ls(filename)
  # extract objects as writtten by h5_write_misc
  if("/genesets" %in% content$group){
    genesets <- list()
    for(name in content$name[which(content$group == "/genesets")]){
      genesets[[name]] <- rhdf5::h5read(filename, paste("genesets/",name,sep="/"))
    }
  }else{
    genesets <- NULL
  }
  algorithms <- rhdf5::h5read(filename, "algorithms")
  grouping <- as.factor(rhdf5::h5read(filename, "grouping"))
  function.call <- as.list(rhdf5::h5read(filename, "function_call/args"))
  names(function.call) <- rhdf5::h5read(filename, "function_call/argnames")
  function.call <- as.call(function.call)
  return(list(genesets = genesets, algorithms = algorithms, grouping = grouping, function.call = function.call))
}
