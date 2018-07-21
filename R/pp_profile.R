##' Plot profiles
##'
##' Plot a phylogenetic profile, which is a matrix.
##'
##' @title pp plot profiles
##' @param x A numeric matrix.
##' @inheritParams pp_tile
##' @return A \code{gg} class object
##' @examples
##' require('ggplot2')
##' require('magrittr')
##'
##' ## a binning profile
##' ppB <- sample(0:1, 5 * 20, replace = TRUE) %>% matrix(nrow = 5)
##' pp_profile(ppB) + scale_fill_manual(values = c('grey91', 'steelblue'))
##'
##' ## a continuous profile
##' ppC <- rnorm(5 * 20) %>% matrix(nrow = 5)
##' pp_profile(ppC) + scale_fill_gradient2()
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom ggplot2 ggplot geom_tile labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @seealso \code{\link[ggplot2]{geom_tile}}
##' @seealso \code{\link[ggplot2]{theme}}
##' @export
##' 
pp_profile <- function(x, legend.position = 'none', ...) {

  nr <- nrow(x)
  nc <- ncol(x)
  if (isBinMat_internal(x)) {
    x <- factor(c(x))
  } else {
    x <- c(x)
  }
  m <- data.frame(x = rep(1:nc, each = nr),
                  y = rep(1:nr, nc),
                  fill = x)

  pObj <- ggplot(m, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fill'), ...) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    theme_pp(legend.position = legend.position)

  return(pObj)
}
