##' Plot tiles
##'
##' Plot tiles to discrete groups
##'
##' \code{pp_htext()}: Plot tiles in the horizontal direction.
##'
##' \code{pp_vtext()}: Plot tiles in the vertical direction.
##'
##' @title pp plot tiles
##' @param x A character vector
##' @param ... Parameters passed to \code{geom_tile()} in the ggplot2 package.
##' @return A \code{gg} class object
##' @examples
##' pp_htile(rep(letters[1:3], 4))
##' pp_vtile(rep(LETTERS[1:2], 3)) + scale_fill_manual(values = c('red','blue', 'green'))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_tile labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @rdname pp_tile
##' @seealso \code{\link[ggplot2]{geom_tile}}
##' @export
##' 
pp_htile <- function(x, ...) {

  xlen <- length(x)
  m <- data.frame(x = rep(0, xlen),
                  y = 1:xlen,
                  fill = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fill'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    theme_pp(legend.position='none')

  return(tObj)
}


##' @inheritParams pp_vtile
##' @rdname pp_tile
##' @export
##' 
pp_vtile <- function(x, ...) {

  xlen <- length(x)
  m <- data.frame(y = rep(0, xlen),
                  x = 1:xlen,
                  fill = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fill'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    theme_pp(legend.position='none')

  return(tObj)
}
