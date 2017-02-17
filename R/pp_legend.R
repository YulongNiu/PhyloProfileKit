##' ggplot legend
##'
##' Retrieve the ggplot legend from a given ggplot2.
##'
##' @title Retrieve ggplot2 legend
##' @param g \code{ggplot2} object.
##' @param ... Parameters passed to \code{theme()} in the ggplot2 package.
##' @return A \code{gtable} object.
##' @examples
##' require('ggplot')
##'
##' p <- ggplot(mtcars, aes(mpg, wt)) +
##'   geom_point(aes(colour = factor(cyl)))
##' plot(pp_legend(p, legend.position = 'left'))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplotGrob theme
##' @importFrom grid unit
##' @references \url{https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs}
##' @references \url{https://stackoverflow.com/questions/17012518/why-does-this-r-ggplot2-code-bring-up-a-blank-display-device}
##' @references \url{https://github.com/hadley/ggplot2/issues/809}
##' @export
##'
pp_legend <- function(g, ...) {
  g <- ggplotGrob(g + theme(legend.spacing = unit(0, 'mm'), ...))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

  return(legend)
}

