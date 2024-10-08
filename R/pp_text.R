##' Plot text
##'
##' Plot text, for example column names and row names.
##'
##' @title pp plot text
##' @param x A character vector
##' @param shift A numeric value indicating shit scale.
##' @param ... Parameters passed to \code{geom_text()} in the ggplot2 package.
##' @inheritParams pp_tile
##' @return A \code{gg} class object
##' @examples
##' require(ggplot2)
##'
##' pp_text(letters[1:10], colour = 'grey55', size = 3)
##' pp_text(rep(LETTERS[1:2], 3), colour = factor(rep(c('blue', 'red'), 3)))
##' pp_text(letters[1:5]) + coord_flip()
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom ggplot2 ggplot geom_text labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @seealso \code{\link[ggplot2]{geom_text}}
##' @seealso \code{\link[ggplot2]{theme}}
##' @export
##' 
pp_text <- function(x, shift = 0.5, legend.position = 'none', ...) {


  xlen <- length(x)
  m <- data.frame(x = rep(0, xlen),
                  y = seq(shift, xlen - shift, 1),
                  label = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_text(aes_string(label = 'label'), ...) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, xlen), breaks = NULL) +
    theme_pp(legend.position = legend.position)

  return(tObj)
}

