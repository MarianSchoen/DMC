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
  content <- h5ls(filename)
  # extract objects as writtten by h5_write_misc
  if("/genesets" %in% content$group){
    genesets <- list()
    for(name in content$name[which(content$group == "/genesets")]){
      genesets[[name]] <- h5read(filename, paste("genesets/",name,sep="/"))
    }
  }else{
    genesets <- NULL
  }
  algorithms <- h5read(filename, "algorithms")
  grouping <- as.factor(h5read(filename, "grouping"))
  function.call <- as.list(h5read(filename, "function_call/args"))
  names(function.call) <- h5read(filename, "function_call/argnames")
  function.call <- as.call(function.call)
  return(list(genesets = genesets, algorithms = algorithms, grouping = grouping, function.call = function.call))
}