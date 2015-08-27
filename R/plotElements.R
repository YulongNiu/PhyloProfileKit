##' Phylo theme
##'
##' A totally blank gpplot2 theme expect for legend
##' @title Blank theme
##' @param ... ggplot2 theme() additional parameters
##' @return ggplot2 object
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 element_blank theme_bw %+replace% theme
##' @importFrom grid unit
##' @keywords internal
##' 
theme_phylo <- function(...) {
  theme_bw() %+replace%
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, 'mm'),
        axis.ticks.margin = unit(0, 'mm'),
        axis.line = element_blank(),
        panel.margin = unit(0, 'mm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'line'),
        legend.margin = unit(0, 'mm'),
        ...)
}

##' Empty white block
##'
##' A empty white block used as interval to organize plot
##' @title Empty block
##' @param ... ggplot2 geom_point() additional parameters
##' @return ggplot2 object
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_point labs aes_string scale_y_continuous scale_x_continuous
##' @keywords internal
##' 
EmptyEle <- function(...) {
  emptyData <- data.frame(x = 1, y = 1)
  empty <- ggplot(emptyData) +
    geom_point(aes_string('x', 'y'), colour='white', ...) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            theme_phylo(legend.position='none')

  return(empty)
}
