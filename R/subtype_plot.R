subtype_plot <- function(results.df){
  
  results.df <- results.df[which(results.df$cell_type == "overall"),]
  p <- ggplot(results.df, aes(x = cluster_size, y = score, col = algorithm)) + 
    geom_line(aes(linetype = coarse)) +
    scale_x_reverse()
  
  return(p)
}