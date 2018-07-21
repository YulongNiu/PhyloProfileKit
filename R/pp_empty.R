##' Plot empty block
##'
##' Plot empty blocks, for example an empty white block.
##'
##' @title pp plot empty
##' @param ... Parameters passed to \code{geom_point()} in the ggplot2 package.
##' @return A \code{gg} class object
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom ggplot2 ggplot geom_point labs aes_string scale_y_continuous scale_x_continuous
##' @examples
##' require(ggplot2)
##'
##' pp_empty(colour = 'white')
##' @export
##' 
pp_empty<- function(...) {
  m <- data.frame(x = 1, y = 1)
  eObj <- ggplot(m) +
    geom_point(aes_string('x', 'y'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    theme_pp(legend.position='none')

  return(eObj)
}
