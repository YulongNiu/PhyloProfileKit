
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
`+.gtable` <- function(x, y) {
  return(cbind(x, y))
}


##' Add plot object to \code{pmat}.
##'
##' %@<%: at left.
##'
##' %@>%: at right.
##'
##' %@^%: at top.
##'
##' %@v%: at bottom.
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
`%@<%` <- function(x, y) {
  x@left <- y
  return(x)
}


##' @inheritParams "%@<%"
##' @rdname loc
##' @export
##' 
`%@>%` <- function(x, y) {
  x@right <- y
  return(x)
}

##' @inheritParams "%@<%"
##' @rdname loc
##' @export
##' 
`%@^%` <- function(x, y) {
  x@top <- y
  return(x)
}

##' @inheritParams "%@<%"
##' @rdname loc
##' @export
`%@v%` <- function(x, y) {
  x@bottom <- y
  return(x)
}


