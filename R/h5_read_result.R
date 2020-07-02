#' read_result_list
#'
#' @param filename string, which file should be read in?
#' @param content h5 file content; used instead of filename in case of recursive function call
#' @param groupname name of the current group
#' @return list
#' @export 
read_result_list <- function(filename, content = NULL, groupname = NULL) {
  # parameter check
  if(!file.exists(filename)) {
    stop(paste("Could not find file ", filename, sep = ""))
  }
  # if content is NULL, read from file
  if(is.null(content)){
    content <- rhdf5::h5dump(filename)
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
    # reconstruct the objects stored in the lowest levels of the list
    if("data" %in% names(content)){
      temp <- as.matrix(content[["data"]])
      # assign dimnames according to the data (groupname)
      if(!is.null(groupname)){
        if(groupname == "est.props" || groupname == "bulk.props" || groupname == "c.true" || groupname == "c.true.coarsly"){
          rownames(temp) <- content[["celltypeids"]]
          colnames(temp) <- content[["bulkids"]]
        }else{
          if(groupname == "sig.matrix" || groupname == "ref.profiles"){
            rownames(temp) <- content[["geneids"]]
            colnames(temp) <- content[["celltypeids"]]
          }else{
            # fallback
            rownames(temp) <- content[["celltypeids"]]
            colnames(temp) <- content[["bulkids"]]
          }
        }
      }
      return(temp)
    }else{
      if(!is.null(groupname) && groupname == "g"){
        temp <- content[["values"]]
        names(temp) <- content[["geneids"]]
        return(temp)
      }
    }
    # if there is no est.props or sig.matrix (lists) but
    # name and times, then est.props and sig.matrix are NULL
    if(all(c("name", "times") %in% names(content))){
      content["est.props"] <- list(NULL)
      content["sig.matrix"] <- list(NULL)
    }
  }
  return(content)
}
