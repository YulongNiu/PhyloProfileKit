
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
`+.gtable` <- function(x, y) {
  return(cbind(x, y))
}


## ##' Add plot object to \code{gmat}.
## ##'
## ##' \code{atleft()}: at left.
## ##'
## ##' \code{atright()}: at right.
## ##'
## ##' \code{attop()}: at top.
## ##'
## ##' \code{atbottom()}: at bottom.
## ##' 
## ##' @title Add locations
## ##' @inheritParams atleft
## ##' @return A \code{gmat} object.
## ##' @examples
## ##' require('gridExtra')
## ##' require('ggplot2')
## ##'
## ##' p <- qplot(1,1)
## ##' new('gmat') %@<% grid.arrange(p)
## ##' @author Yulong Niu \email{niuylscu@@gmail.com}
## ##' @rdname location-methods
## ##' @exportMethod atleft
## ##'
## setMethod(f = 'atleft',
##           signature = c(x = 'gmat', y = 'gtable'),
##           definition = function(x, y, ...) {
##             x@left <- y
##             return(x)
##           })

## ##' @inheritParams atleft
## ##' @rdname location-methods
## ##' @exportMethod atright
## ##'
## setMethod(f = 'atright',
##           signature = c(x = 'gmat', y = 'gtable'),
##           definition = function(x, y, ...) {
##             x@right <- y
##             return(x)
##           })

## ##' @inheritParams atleft
## ##' @rdname location-methods
## ##' @exportMethod attop
## ##'
## setMethod(f = 'attop',
##           signature = c(x = 'gmat', y = 'gtable'),
##           definition = function(x, y, ...) {
##             x@top <- y
##             return(x)
##           })

## ##' @inheritParams atleft
## ##' @rdname location-methods
## ##' @exportMethod atbottom
## ##' 
## setMethod(f = 'atbottom',
##           signature = c(x = 'gmat', y = 'gtable'),
##           definition = function(x, y, ...) {
##             x@bottom <- y
##             return(x)
##           })


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

