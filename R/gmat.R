##' Add \code{ggplot} and \code{gtable} objects
##'
##' @title Add methods
##' @param e1 A \code{ggplot} or code{gtable} object.
##' @param e2 A \code{ggplot} or code{gtable} object.
##' @return A \code{gtable} object
##' @examples
##' require('ggplot2')
##' require('gridExtra')
##'
##' p <- ggplot(mtcars, aes(mpg, wt)) +
##'   geom_point(aes(colour = factor(cyl)))
##' p %@+% p
##' p %@+% grid.arrange(qplot(1,1)) %@+% p
##' p %@+% (p + scale_colour_manual(values = c('red', 'blue', 'green')))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom gridExtra grid.arrange
##' @rdname addinloc
##' @export
##' 
atadd <- function(e1, e2) {

  ## check obj is ggplot or gtable
  checkg_internal <- function(x) {
    ifelse(inherits(x, 'ggplot') ||
           inherits(x, 'gtable'),
           TRUE,
           FALSE)
  }

  stopifnot(checkg_internal(e1) &&
            checkg_internal(e2))

  e1 <- gg2t(e1)
  e2 <- gg2t(e2)

  return(cbind(e1, e2))
}

##' @rdname addinloc
##' @export
##'
`%@+%` <- atadd

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
##' @param y A \code{ptable} or \code{ggplot} object.
##' @return A \code{pmat} object.
##' @examples
##' require('gridExtra')
##' require('ggplot2')
##'
##' p <- qplot(1,1)
##' new('gmat') %@<% p
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname loc
##' @export
##'
atleft <- function(x, y) {
  x@left <- gg2t(y)
  return(x)
}


##' @inheritParams atleft
##' @rdname loc
##' @export
##'
atright <- function(x, y) {
  x@right <- gg2t(y)
  return(x)
}

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
attop <- function(x, y) {
  x@top <- gg2t(y)
  return(x)
}

##' @inheritParams atleft
##' @rdname loc
##' @export
atbottom <- function(x, y) {
  x@bottom <- gg2t(y)
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


##' Transfer to \code{gtable}
##'
##' @title Transfer to \code{gtable} objects.
##' @param x A \code{ggplot} or \code{gtable} object.
##' @return A \code{gtable} object.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom gridExtra grid.arrange
##' @keywords internal
##' 
gg2t <- function(x) {
  if (inherits(x, 'ggplot')) {
    x <- grid.arrange(x)
  } else {}

  return(x)
}



