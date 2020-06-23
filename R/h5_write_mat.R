#' write a simple matrix to a hdf5 file
#' 
#' @param m matrix with > 0 rows and > 0 named columns
#' @param filename filename of hdf5 file to be written to
#' @param filepath optional, path where the file should be stored.
#' Only necessary, if the directory does not yet exist, otherwise it can
#' be supplied as part of filename
#' @return NULL, write to file

h5_write_mat <- function(m, filename, filepath = ""){
    # check input parameters
    if(is.null(m) || nrow(m) == 0 || ncol(m) == 0){
        warning("Invalid data frame. Cannot write to file.")
        return(NULL)
    }
    if(length(colnames(m)) != ncol(m)){
        warning("matrix needs to have proper column names.")
        return(NULL)
    }
    if(!is.character(filename) || !is.character(filepath)){
        warning("Invalid file path. Cannot write to file.")
        return(NULL)
    }
    if(filepath != ""){
        if(!dir.exists(filepath)){
            flag <- dir.create(filepath, recursive = TRUE)
            if(!flag){
                warning("Could not create file path. Return without writing.")
                return(NULL)
            }
        }
    }

    # create h5 file
    if(!file.exists(filename)){
	    rhdf5::h5createFile(filename)
    }

    # write row names
    if(!is.null(rownames(m)) && !is.na(rownames(m))){
	    rhdf5::h5write(as.vector(rownames(m)), filename, "rownames", write.attributes = TRUE)
    }
    # write column names
    rhdf5::h5write(as.vector(colnames(m)), filename, "colnames", write.attributes = TRUE)
    # write matrix
    rhdf5::h5write(as.matrix(m), filename, "data", write.attributes = FALSE)
    
}
