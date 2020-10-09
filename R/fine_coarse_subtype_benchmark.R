#' fine and coarse subtype benchmark function
#'
#' # basically, there are 3 different things to compare. 
#' 1. correlation per algorithm on the deep C 
#' 2. correlation per algorithm on the accumulated cursory C 
#' 3. special scenario where subtype 1 is in X, and subtype 2 is in Y
#'
#' 1. and 2. can be done in one go, 3. needs more prepration. 
#'
#' @param sc.counts 
#' @param sc.pheno 
#' @param subtype.pattern 
#' @param cell.type.column 
#' @param sample.name.column 
#' @param avg.profiles.per.subcluster 
#' @param verbose 
#' @param algorithm.list 
#'
#' @return
#'
#' @examples
fine_coarse_subtype_benchmark <- function(
  sc.counts = cll.exprs, 
  sc.pheno = cll.pheno, 
  subtype.pattern = "subtype",
  cell.type.column = "cell_type", 
  sample.name.column = "sample.name", 
  avg.profiles.per.subcluster =  NULL,
  n.cluster.sizes = 5,
  verbose = TRUE, 
  algorithm.list = list(
    list(algorithm = run_dtd, name = "DTD"),
    list(algorithm = run_least_squares, name = "Least_Squares")
  ),
  patient.column = NULL,
  n.bulks = 500
){
  if(is.null(avg.profiles.per.subcluster)){
    min.profiles <- max(2, as.integer(min(table(sc.pheno[[cell.type.column]])) / 2) - 1)
    max.profiles <- as.integer(max(table(sc.pheno[[cell.type.column]])) / 2) + 1

    avg.profiles.per.subcluster <- as.integer(seq(min.profiles, max.profiles, length.out = 5))
  }else{
    n.cluster.sizes <- length(avg.profiles.per.subcluster)
  }
  if(ncol(sc.counts) != nrow(sc.pheno)){
    stop(
      "In DAB::fine_coarse_subtype_benchmark: 'sc.counts' does not fit 'sc.pheno'
      'ncol(sc.counts)' must equal 'nrow(sc.pheno)'"
    )
  }
  
  if(! cell.type.column %in% colnames(sc.pheno)){
    stop(
      "In DAB::fine_coarse_subtype_benchmark: 'cell.type.column' is not in 
      'colnames(sc.pheno)'"
    )
  }  
  
  if(! sample.name.column %in% colnames(sc.pheno)){
    stop(
      "In DAB::fine_coarse_subtype_benchmark: 'sample.name.column' is not in 
      'colnames(sc.pheno)'"
    )
  }  
  
  # find a colname for the subtype information:   
  while (
    any(
      grepl(
        pattern = subtype.pattern
        , x = colnames(sc.pheno)
      )
    )
  ){
    subtype.pattern <- paste0(subtype.pattern, sample.int(100, 1))
  }
  
  # subcluster the data, return the sc.pheno frame with additional columns
  sc.pheno <- add_subtype_pheno(
    sc.counts = sc.counts
    , sc.pheno = sc.pheno
    , cell.type.column = cell.type.column
    , sample.name.column = sample.name.column
    , new.subtype.column = subtype.pattern
    , hclust.obj = NA
    , avg.profiles.per.subcluster.vec = avg.profiles.per.subcluster
    , verbose = verbose
  )
  
  # detect the "new" subtype columns:
  subtype.columns <- colnames(sc.pheno)[
    grepl(
      pattern = subtype.pattern
      , x = colnames(sc.pheno)
    )
    ]
  
 
  # initialise the output list. 
  # (this holds as many entries as there are 'subtype.columns')
  # each entry holds the following entries: 
  #  1. a list entry for each 'subtype.column'
  #  2. There are four entries: 
  #   - 'c.true': matrix, true cell compositions, for all subtypes
  #   - 'c.true.coarsly': matrix, true cell compositions, for the 
  #                       'cell.type.column'cell types
  #   - 'c.estimated.list': a list, with a matrix entry for each algorithm. 
  #          There, the estimated cell compositions, for all subtypes are stored
  #   - 'c.estimated.coarsly.list': a list, with a matrix entry for each algorithm. 
  #          There, the estimated cell compositions, for the 
  #          'cell.type.column'cell types are stored.
  column.list <- vector(
    mode = "list"
    , length = length(subtype.columns)
  )
  names(column.list) <- subtype.columns
 
  # add original cell type classification as new column for benchmark
  #profiles.per.ct <- as.integer(max(table(sc.pheno[[cell.type.column]]))/2) + 1
	  #as.integer(nrow(sc.pheno) / length(unique(sc.pheno[[cell.type.column]]))) 
  #sc.pheno <- cbind(sc.pheno, sc.pheno[[cell.type.column]])
  #colnames(sc.pheno)[ncol(sc.pheno)] <- paste0("subtype.avg.", profiles.per.ct)
  #subtype.columns <- c(paste0("subtype.avg.", profiles.per.ct), subtype.columns)

  # go through each 'subtype.column' 
  for(column in subtype.columns){
    # simulate bulks, and deconvolute them 
    # (currently, I don't seperate test/train)
    if(verbose) cat(column, "\n")
    some.estimates <- deconvolute(
      training.expr = sc.counts
      , training.pheno = sc.pheno
      , test.expr = sc.counts
      , test.pheno = sc.pheno
      , cell.type.column = column
      , algorithms = algorithm.list
      , verbose = verbose
      , split.data = FALSE
      , n.repeats = 1 # don't yet increase this please
      , patient.column = patient.column
      , n.bulks = n.bulks
    ) 
    
    # extract the true c (with all subtype quantities)... 
    c.true <- some.estimates$bulk.props
    # this C is quite deep. I need the coarse C, therefore:
    c.true.coarsly <- matrix(
      data = NA
      , nrow = length(unique(sc.pheno[, cell.type.column]))
      , ncol = ncol(c.true)
    )
    colnames(c.true.coarsly) <- colnames(c.true)
    rownames(c.true.coarsly) <- unique(sc.pheno[, cell.type.column])
    for(major.cell.type in rownames(c.true.coarsly)){
      associated.subtypes <- rownames(c.true)[
        which(
          grepl(
            pattern = major.cell.type
            , x = rownames(c.true)
          )
        )
        ]
        if(length(associated.subtypes) > 0){
          c.true.coarsly[major.cell.type, ] <- colSums(
                  x = c.true[associated.subtypes, , drop = FALSE]
                )
        }else{
          c.true.coarsly[major.cell.type, ] <- 0
        }
      
    }
    
    # store the true C matrices: 
    column.list[[column]][["c.true"]] <- c.true
    column.list[[column]][["c.true.coarsly"]] <- c.true.coarsly
    
    # TODO: manage the "n.repeats = 1" problem
    algorithms <- names(some.estimates$results.list$`1`)
    
    # for the estimated C, initialise two list entries: 
    column.list[[column]][["c.estimated.list"]] <- vector(
      mode = "list"
      , length = length(algorithms)
    )
    names(column.list[[column]][["c.estimated.list"]]) <- algorithms
    
    
    column.list[[column]][["c.estimated.coarsly.list"]] <- vector(
      mode = "list"
      , length = length(algorithms)
    )
    names(column.list[[column]][["c.estimated.coarsly.list"]]) <- algorithms
    
    # go through all provided algorithms
    for(algorithm in algorithms){
      # for the current algorithm, extract the estimated 'C's
      c.estimated <- some.estimates$results.list$`1`[[algorithm]]$est.props
      # store it: 
      column.list[[column]][["c.estimated.list"]][[algorithm]] <- c.estimated
      
      # add up those subtypes, that origin from the same cell type:       
      c.estimated.coarsly <- matrix(
        data = NA
        , nrow = nrow(c.true.coarsly)
        , ncol = ncol(c.true.coarsly)
      )
      dimnames(c.estimated.coarsly) <- dimnames(c.true.coarsly)
      for(major.cell.type in rownames(c.true.coarsly)){
        associated.subtypes <- rownames(c.estimated)[
          which(
            grepl(
              pattern = major.cell.type
              , x = rownames(c.estimated)
            )
          )
          ]
          if(length(associated.subtypes) > 0){
            c.estimated.coarsly[major.cell.type, ] <- colSums(
              x = c.estimated[associated.subtypes, , drop = FALSE]
            )
          }else{
            c.estimated.coarsly[major.cell.type, ] <- 0
          }
        
      }
      # and store again
      column.list[[column]][["c.estimated.coarsly.list"]][[algorithm]] <- c.estimated.coarsly
    }
  }
  return(column.list)
}
