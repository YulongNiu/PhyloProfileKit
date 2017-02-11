
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
`+.gtable` <- function(x, y) {
  return(cbind(x, y))
}

##' Add plot object to \code{pmat}.
##'
##' \code{atleft()}/%@<%: at left.
##'
##' \code{atright()}/%@>%: at right.
##'
##' \code{attop()}/%@^%: at top.
##'
##' \code{atbottom()}/%@v%: at bottom.
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

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
`%@<%` <- atleft

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
`%@>%` <- atright

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
`%@^%` <- attop

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
`%@v%` <- atbottom



