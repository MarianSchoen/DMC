#' read a simple matrix from hdf5 file
#' 
#' @param filename character specifying the file to read from
#' @return matrix read from file

h5_read_mat <- function(filename){
    # check input parameter
    if(!is.character(filename) || !file.exists(filename)){
        warning("Invalid filename, returning NULL")
        return(NULL)
    }

    # list the content and load the data
    content <- h5ls(filename)
    if("rownames" %in%  content$name){
        rnames <- h5read(filename, "rownames")
    }
    cnames <- h5read(filename, "colnames")
    m <- h5read(filename, "data")
    colnames(m) <- cnames
    if(exists("rnames")){
        rownames(m) <- rnames
    }
    
    return(m)
}