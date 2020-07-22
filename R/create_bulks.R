#' simulate bulk data by summing over single-cell data
#'
#' @param exprs non-negative numeric, scRNA-seq profiles as columns, 
#' features as rows
#' @param pheno data.frame, phenotype data containing
#' cell type labels for the expression matrix,
#' must contain 'cell.type.column',
#' ordering of cells must be the same as in exprs
#' @param n.bulks integer, the number of bulks to be created, defaults to 500
#' @param include.in.bulks vector of strings, cell types to be used for bulk 
#' simulation; if not supplied, all will be used
#' @param fraction.per.bulk positive numeric < 1, fraction of samples to be randomly
#' drawn for each bulk; default 0.1
#' @param sum.to.count boolean, should all bulks be normalized
#' to a fixed total count number? default TRUE
#' @param cell.type.column string, which column of 'pheno'
#' holds the cell type information? 
#' @return list with 
#'    - "bulks": matrix containing bulk expression profiles (features x bulks)
#'    - "props" matrix containing quantities (cell type x bulks)
#'    - "sub.props": TODO
#' @example create_bulks(training.exprs, training.pheno, n.bulks = 1000)
create_bulks <- function(
    exprs, 
    pheno, 
    cell.type.column = "cell_type",
    n.bulks = 500, 
    include.in.bulks = NULL, 
    fraction.per.bulk = 0.1, 
    sum.to.count = TRUE
    ) {
    # parameter checks
    if (nrow(pheno) != ncol(exprs)) {
        stop("Number of columns in exprs and rows in pheno do not match")
    }
    if(!is.numeric(n.bulks)){
        stop("n.bulks must be numeric")
    }
    if(!is.numeric(fraction.per.bulk)){
        stop("fraction.per.bulk must be numeric (0,1)")
    }else{
        if(fraction.per.bulk <= 0 || fraction.per.bulk >= 1){
            stop("fraction.per.bulk must be numeric (0,1)")
        }
    }
    rownames(pheno) <- colnames(exprs)

    # keep only specified cell types
    if (is.null(include.in.bulks)) {
        include.in.bulks <- unique(pheno[, cell.type.column])
    }

    if (length(which(pheno[, cell.type.column] %in% include.in.bulks)) > 0) {
        pheno <- pheno[which(pheno[, cell.type.column] %in% include.in.bulks), , drop = FALSE]
    }
    exprs <- exprs[, rownames(pheno), drop = FALSE]

    # scale to fixed total count per profile
    exprs <- scale_to_count(exprs)

    # create a matrix to contain the bulk expression profiles
    bulk.exprs <- matrix(
        0,
        nrow = nrow(exprs),
        ncol = n.bulks
    )
    rownames(bulk.exprs) <- rownames(exprs)
    colnames(bulk.exprs) <- as.character(1:n.bulks)

    # create a matrix to contain true proportions for each simulated bulk
    props <- matrix(0, nrow = length(unique(pheno[, cell.type.column])), ncol = n.bulks)
    rownames(props) <- unique(pheno[, cell.type.column])
    colnames(props) <- colnames(bulk.exprs)
    
    # if column subtypes exists create a matrix containing proportions of subtypes
    if("subtype" %in% colnames(pheno)){
        subtypes <- unique(pheno$subtype)
        n.subtypes <- length(subtypes)
        combined.type <- paste(pheno[[cell.type.column]], pheno$subtype, sep = ".")
        sub.props <- matrix(0, length(unique(combined.type)), ncol = n.bulks)
        rownames(sub.props) <- unique(combined.type)
        colnames(sub.props) <- colnames(bulk.exprs)
    }else{
        sub.props <- NULL
        n.subtypes <- 0
    }

    # if no subtype, sample with random distribution from cell_types
    # if subtypes, sample with random distribution from combined.type
    
    for(i in 1:ncol(bulk.exprs)){
      bulk.samples <- c()
      if("subtype" %in% colnames(pheno)){
        types <- unique(combined.type)
        
        # sample random weights for each type
        weights <- sample(0:100, size = length(types), replace = TRUE)
        
        # determine, how many samples of each type should be drawn
        n.samples <- sample(types, size = ceiling(fraction.per.bulk * ncol(exprs)), replace = TRUE, prob = weights)
        
        # draw samples for each type, store the proportions and the expression
        for(t in types){
          sub.samples <- sample(which(combined.type == t), size = sum(n.samples == t), replace = TRUE)
          coarse.types <- pheno[[cell.type.column]][sub.samples]
          bulk.samples <- c(bulk.samples, sub.samples)
          
          sub.props[t, i] <- length(sub.samples)
          for(ct in unique(coarse.types)){
            props[ct, i] <- props[ct, i] + sum(coarse.types == ct)
          }
        }
        sub.props[,i] <- sub.props[,i] / length(bulk.samples)
        props[,i] <- props[,i] / length(bulk.samples)
      }else{
        types <- unique(pheno[, cell.type.column])
        # sample random weights for each type
        weights <- sample(0:100, size = length(types), replace = TRUE) #
        # determine, how many samples of each type should be drawn
        n.samples <- sample(types, size = ceiling(fraction.per.bulk * ncol(exprs)), replace = TRUE, prob = weights) #
        
        # draw samples for each type, store the proportions and the expression
        for(t in types){
          coarse.samples <- sample(which(pheno[,cell.type.column] == t), size = sum(n.samples == t), replace = TRUE)
          bulk.samples <- c(bulk.samples, coarse.samples)
          props[t, i] <- length(coarse.samples)
        }
        props[,i] <- props[,i] / length(bulk.samples)
      }
      
      bulk.exprs[,i] <- rowSums(exprs[, bulk.samples, drop = F])
    }

    # let no feature occur twice in the bulks
    if(any(duplicated(rownames(bulk.exprs)))){
        warning("Found duplicate features in simulated bulks. Removing...")
        bulk.exprs <- bulk.exprs[-which(duplicated(bulk.exprs)), ]
    }

    # sum bulks to fixed total count per profile if sum.to.count is true
    if (sum.to.count) {
        bulk.exprs <- apply(bulk.exprs, 2, function(x) {
            (x / sum(x)) * length(x)
        })
    }

    return(list(bulks = bulk.exprs, props = props, sub.props = sub.props))
}
