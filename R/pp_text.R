##' Plot text
##'
##' Plot text, for example column names and row names.
##'
##' \code{pp_htext()}: Plot text in the horizontal direction.
##'
##' \code{pp_vtext()}: Plot text in the vertical direction.
##'
##' @title pp plot text
##' @param x A character vector
##' @param shift A numeric value indicating shit scale.
##' @param ... Parameters passed to \code{geom_text()} in the ggplot2 package.
##' @return A \code{gg} class object
##' @examples
##' pp_htext(letters[1:10], colour = 'grey55', size = 3)
##' pp_vtext(rep(LETTERS[1:2], 3), colour = factor(rep(c('blue', 'red'), 3)))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_text labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @rdname plottext
##' @seealso \code{\link[ggplot2]{geom_text}}
##' @export
##' 
pp_htext <- function(x, shift = 0.5, ...) {

  xlen <- length(x)
  m <- data.frame(x = rep(0, xlen),
                  y = seq(shift, xlen - shift, 1),
                  label = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_text(aes_string(label = 'label'), ...) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, xlen), breaks = NULL) +
    theme_pp(legend.position='none')

  return(tObj)
}


##' @inheritParams pp_htext
##' @rdname plottext
##' @export
##' 
pp_vtext <- function(x, shift = 0.5, ...) {
  xlen <- length(x)
  m <- data.frame(y = rep(0, xlen),
                  x = seq(shift, xlen - shift, 1),
                  label = x)

  tObj <- ggplot(m, aes_string('x', 'y')) +
    geom_text(aes_string(label = 'label'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, xlen), breaks = NULL) +
    theme_pp(legend.position='none')

  return(tObj)
}

