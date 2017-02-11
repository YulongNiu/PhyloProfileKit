##' Add \code{ggplot} and \code{gtable} objects
##'
##' @title Add methods
##' @param e1 A \code{ggplot} or code{gtable} object.
##' @param e2 A \code{ggplot} or code{gtable} object.
##' @return A \code{gtable} object
##' @examples
##' require('ggplot2')
##'
##' p <- ggplot(mtcars, aes(mpg, wt)) +
##'   geom_point(aes(colour = factor(cyl)))
##' p + p
##' p + (p + scale_colour_manual(values = c("red","blue", "green")))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom gridExtra grid.arrange
##' @importFrom gglot2 ggplotGrob
##' @rdname add-methods
##' @exportMethod '+'
##' 
setMethod(f = '+',
          signature = c(e1 = 'ANY', e2 = 'ANY'),
          definition = function(e1, e2) {

            ## check obj is ggplot or gtable
            checkg_internal <- function(x) {
              ifelse(inherits(x, 'ggplot') ||
                     inherits(x, 'gtable'),
                     TRUE,
                     FALSE)
            }

            stopifnot(checkg_internal(e1) &&
                      checkg_internal(e2))

            if (inherits(e1, 'ggplot')) {
              e1 <- ggplotGrob(e1)
            }
            else if (inherits(e2, 'ggplot')) {
              e2 <- ggplotGrob(e2)
            }
            else {}

            return(grid.arrange(e1, e2, ncol = 2))
          })

##' Add plot object to \code{pmat}.
##'
##' \code{atleft()}: at left.
##'
##' \code{atright()}: at right.
##'
##' \code{attop()}: at top.
##'
##' \code{atbottom()}: at bottom.
##'
##' @title Add locations
##' @param x A \code{pmat} object.
##' @param y A \code{ptable} object.
##' @return A \code{pmat} object.
##' @examples
##' require('gridExtra')
##' require('ggplot2')
##'
##' p <- qplot(1,1)
##' new('gmat') %@<% grid.arrange(p)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname loc
##' @export
##'
atleft <- function(x, y) {
  x@left <- y
  return(x)
}


##' @inheritParams atleft
##' @rdname loc
##' @export
##'
atright <- function(x, y) {
  x@right <- y
  return(x)
}

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
attop <- function(x, y) {
  x@top <- y
  return(x)
}

##' @inheritParams atleft
##' @rdname loc
##' @export
atbottom <- function(x, y) {
  x@bottom <- y
  return(x)
}

##' @rdname loc
##' @export
##'
`%@<%` <- atleft

##' @rdname loc
##' @export
##'
`%@>%` <- atright

##' @rdname loc
##' @export
##'
`%@^%` <- attop

##' @rdname loc
##' @export
##'
`%@v%` <- atbottom

