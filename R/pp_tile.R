##' Plot tiles
##'
##' Plot tiles to discrete groups
##'
##' @title pp plot tiles
##' @param x A character vector
##' @param legend.position Position of the legend and the default is "none" (no legend). See the \code{theme()} in the ggplot2 package.
##' @param ... Parameters passed to \code{geom_tile()} in the ggplot2 package.
##' @return A \code{gg} class object
##' @examples
##' require('ggplot2')
##'
##' pp_tile(rep(letters[1:3], 4))
##' pp_tile(rep(LETTERS[1:2], 3)) + scale_fill_manual(values = c('red','blue', 'green'))
##' pp_tile(rep(letters[1:3], 4)) + coord_flip()
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_tile labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @seealso \code{\link[ggplot2]{geom_tile}}
##' @seealso \code{\link[ggplot2]{theme}}
##' @export
##' 
pp_tile <- function(x, legend.position = 'none', ...) {

  xlen <- length(x)
  m <- data.frame(x = rep(0, xlen),
                  y = 1:xlen,
                  fill = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fill'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    theme_pp(legend.position = legend.position)

  return(tObj)
}
