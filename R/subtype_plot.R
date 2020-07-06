subtype_plot <- function(results.df){
  
  results.df <- results.df[which(results.df$cell_type == "overall"),]
  overall.plot <- ggplot(results.df, aes(x = cluster_size, y = score, col = algorithm)) + 
    geom_line(aes(linetype = coarse)) +
    scale_x_reverse()
  algorithm.plot <- ggplot(results.df, aes(x = cluster_size, y = score, col = algorithm)) + 
    geom_line(aes(linetype = coarse)) +
    scale_x_reverse() +
    facet_grid(rows = vars(algorithm))
  coarse.plot <- ggplot(results.df, aes(x = cluster_size, y = score, col = algorithm)) +
    scale_x_reverse() +
    geom_line() +
    facet_grid(cols = vars(coarse))

  return(list(overall = overall.plot, algorithm = algorithm.plot, coarse = coarse.plot))
}