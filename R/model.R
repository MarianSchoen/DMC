predict_value <- function(df, offset, E){
  result <- 1 / sqrt(1 + abs(E) / df$variance) - abs(offset)
  return(result)
}
