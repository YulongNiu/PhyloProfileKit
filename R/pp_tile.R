pp_tile <- function(x, ...) {

  xlen <- length(x)
  m <- data.frame(x = rep(0, xlen),
                  y = 1:xlen,
                  fill = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fill'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ## scale_fill_manual(values = colour) +
    theme_pp(legend.position='none')

  tGrid <- grid.arrange(tObj, ncol = 1)

  return(tGrid)

}
