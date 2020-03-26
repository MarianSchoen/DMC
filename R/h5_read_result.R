#' read_result_list
#'
#' @param filename string, which file should be read in?
#' @param content 
#' @param groupname 
#'
#' @return list
read_result_list <- function(filename, content = NULL, groupname = NULL) {
  library(rhdf5)
  if(is.null(content)){
    content <- h5dump(filename)
  }
  # go down recursively as in the write function
  # return the content in each function call and add it to the larger list one
  # level above...
  if(any(sapply(content, function(x){is.list(x)}))){
    for(name in names(content)){
      if(is.list(content[[name]])){
        content[[name]] <- read_result_list(filename, content[[name]], name)
      }
    }
  }else{
    if("data" %in% names(content)){
      temp <- as.matrix(content[["data"]])
      if(!is.null(groupname)){
        if(groupname == "est.props" || groupname == "bulk.props"){
          rownames(temp) <- content[["celltypeids"]]
          colnames(temp) <- content[["bulkids"]]
        }else{
          if(groupname == "sig.matrix"){
            rownames(temp) <- content[["geneids"]]
            colnames(temp) <- content[["celltypeids"]]
          }
        }
      }
      return(temp)
    }
    if(all(c("name", "times") %in% names(content))){
      content["est.props"] <- list(NULL)
      content["sig.matrix"] <- list(NULL)
    }
  }
  return(content)
}
