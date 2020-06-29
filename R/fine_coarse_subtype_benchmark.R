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
#' @export
#'
#' @examples
fine_coarse_subtype_benchmark <- function(
  sc.counts = cll.exprs, 
  sc.pheno = cll.pheno, 
  subtype.pattern = "subtype",
  cell.type.column = "cell_type", 
  sample.name.column = "sample.name", 
  avg.profiles.per.subcluster =  c(2, 5, 10, 20, 50, 100), 
  verbose = TRUE, 
  algorithm.list = list(
    list(algorithm = run_dtd, name = "DTD"),
    list(algorithm = run_least_squares, name = "Least_Squares")
  )
){
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
  #  2. then, there are two entries "plot.list", and "cor.list". 
  #  3. there each algorithm has an entry
  #  4a) for the cor.list path there are 2 entries: 
  #      a vector holding per cell type correlations
  #      and a single value averaged over all cell types
  #  4b) for the "plot.list" there is a ggplot object scatter plot for each cell type
  column.list <- vector(
    mode = "list"
    , length = length(subtype.columns)
  )
  names(column.list) <- subtype.columns
  
  # go through each 'subtype.column' 
  for(column in subtype.columns){
    
    # simulate bulks, and deconvolute them 
    # (currently, I don't seperate test/train)
    some.estimates <- deconvolute(
      training.expr = sc.counts
      , training.pheno = sc.pheno
      , test.expr = sc.counts
      , test.pheno = sc.pheno
      , cell.type.column = column
      , algorithms = algorithm.list
      , verbose = TRUE
      , split.data = FALSE
      , n.repeats = 1 # don't yet increase this please
    ) 
    
    # extract the true c
    c.true <- some.estimates$bulk.props
    # this C is quite deep. I need the coarse C, therefore:
    c.true.coarsly <- matrix(
      data = NA
      , nrow = length(unique(sc.pheno[, original.cell.type.column]))
      , ncol = ncol(c.true)
    )
    colnames(c.true.coarsly) <- colnames(c.true)
    rownames(c.true.coarsly) <- unique(sc.pheno[, original.cell.type.column])
    for(major.cell.type in rownames(c.true.coarsly)){
      associated.subtypes <- rownames(c.true)[
        which(
          grepl(
            pattern = major.cell.type
            , x = rownames(c.true)
          )
        )
        ]
      c.true.coarsly[major.cell.type, ] <- colSums(
        x = c.true[associated.subtypes, , drop = FALSE]
      )
    }
    
    # TODO: manage the "n.repeats = 1" problem
    algorithms <- names(some.estimates$results.list$`1`)
    
    # initialise an empty list entry for 'cor.list' and ...
    column.list[[column]][["cor.list"]] <- vector(
      mode = "list"
      , length = length(algorithms)
    )
    names(column.list[[column]][["cor.list"]]) <- algorithms
    # (add an empty vector to 'per.ct'. (because, I c(...) later on.))
    column.list[[column]][["cor.list"]] <- lapply(
      X = column.list[[column]][["cor.list"]]
      , FUN = function(algorithm){
        return(list("per.ct" = c()))
      })
    
    
    
    # ... 'coarse.cor.list', ...
    column.list[[column]][["coarse.cor.list"]] <- column.list[[column]][["cor.list"]]
    
    # ... 'plot.list' ...
    column.list[[column]][["plot.list"]] <- column.list[[column]][["cor.list"]] # false friend, but structure is the same
    
    # and 'coarse.plot.list' ...
    column.list[[column]][["coarse.plot.list"]] <- column.list[[column]][["cor.list"]] # false friend, but structure is the same
    
    
    for(algorithm in algorithms){
      # for the current algorithm, extract the estimated 'C's
      c.estimated <- some.estimates$results.list$`1`[[algorithm]]$est.props
      for(cell.subtype in unique(sc.pheno[, column])){
        a.frame <- data.frame(
          "true" = c.true[cell.subtype, ],
          "estimated" = c.estimated[cell.subtype, ]
        )
        cor.tmp <- cor(a.frame$true,a.frame$estimated)
        column.list[[column]][["cor.list"]][[algorithm]][["per.ct"]][[cell.subtype]] <- cor.tmp
        column.list[[column]][["plot.list"]][[algorithm]][[cell.subtype]] <- ggplot(
          data = a.frame
          , aes(x = true, y = estimated)
        ) + 
          geom_point() + 
          ggtitle(
            paste0(
              cell.subtype, 
              ", Cor: ", 
              round(cor.tmp, digits = 2)
            )
          )
      }
      
      column.list[[column]][["cor.list"]][[algorithm]][["avg.cor"]] <- mean(
        x = column.list[[column]][["cor.list"]][[algorithm]][["per.ct"]]
      )
      
      for(major.cell.type in rownames(c.true.coarsly)){
        associated.subtypes <- rownames(c.estimated)[
          which(
            grepl(
              pattern = major.cell.type
              , x = rownames(c.estimated)
            )
          )
          ]
        a.frame <- data.frame(
          "true" = c.true.coarsly[major.cell.type, ],
          "estimated" = colSums(
            x = c.estimated[associated.subtypes, , drop = FALSE]
          )
        )
        cor.tmp <- cor(a.frame$true, a.frame$estimated)      
        column.list[[column]][["coarse.cor.list"]][[algorithm]][["per.ct"]][[cell.subtype]] <- cor.tmp
        
        column.list[[column]][["coarse.plot.list"]][[algorithm]][[cell.subtype]] <- ggplot(
          data = a.frame
          , aes(x = true, y = estimated)
        ) + 
          geom_point() + 
          ggtitle(
            paste0(
              cell.subtype, 
              ", Cor: ", 
              round(cor.tmp, digits = 2)
            )
          )
      }
      column.list[[column]][["coarse.cor.list"]][[algorithm]][["avg.cor"]] <- mean(
        x = column.list[[column]][["coarse.cor.list"]][[algorithm]][["per.ct"]]
      )
    }
  }
  return(column.list)
}
