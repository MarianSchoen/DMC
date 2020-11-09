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
	if(exists("rnames")){
		rm("rnames")
	}
    # list the content and load the data
    content <- rhdf5::h5ls(filename)
    if("rownames" %in%  content$name){
        rnames <- rhdf5::h5read(filename, "rownames")
    }
    cnames <- rhdf5::h5read(filename, "colnames")
    m <- rhdf5::h5read(filename, "data")
    colnames(m) <- cnames
    if(exists("rnames")){
	    if(length(rnames) == nrow(m)){
        	rownames(m) <- rnames
	    }
    }
    
    return(m)
}
