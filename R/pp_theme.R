##' PhyloProfile theme
##'
##' A totally blank gpplot2 theme expect for the legend
##'
##' @title Blank theme
##' @param ...  Parameters passed to \code{theme()} in the ggplot2 package.
##' @return ggplot2 object
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 element_blank theme_bw %+replace% theme
##' @importFrom grid unit
##' @seealso \code{\link[ggplot2]{theme}}
##' @export
##' 
theme_pp <- function(...) {
  theme_bw() %+replace%
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, 'mm'),
        axis.line = element_blank(),
        panel.spacing = unit(0, 'mm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'line'),
        legend.spacing = unit(0, 'mm'),
        ...)
}
