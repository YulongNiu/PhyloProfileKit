##' @include AllClasses.R AllGenerics.R
NULL

##' Connect \code{ggplot} objects to a list.
##'
##' @title Add methods
##' @param e1 A \code{ggplot} or list.
##' @param e2 A \code{ggplot} or list.
##' @return A list.
##' @examples
##' require('ggplot2')
##' require('gridExtra')
##'
##' p <- ggplot(mtcars, aes(mpg, wt)) +
##'   geom_point(aes(colour = factor(cyl)))
##' p %@+% p
##' q <- p %@+% qplot(1,1) %@+% p
##' p %@+% (p + scale_colour_manual(values = c('red', 'blue', 'green')))
##'
##' ## plot
##' grid.arrange(grobs = q, ncol = 1)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname addinloc
##' @export
##' 
atadd <- function(e1, e2) {

  e1 <- gg2list(e1)
  e2 <- gg2list(e2)

  return(c(e1, e2))
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
##' require('ggplot2')
##'
##' p <- qplot(1,1)
##' new('gmat') %@<% p
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname loc
##' @export
##'
atleft <- function(x, y) {
  x@left <- c(x@left, gg2list(y))
  return(x)
}


##' @inheritParams atleft
##' @rdname loc
##' @export
##'
atright <- function(x, y) {
  x@right <- c(x@right, gg2list(y))
  return(x)
}

##' @inheritParams atleft
##' @rdname loc
##' @export
##'
attop <- function(x, y) {
  x@top <- c(x@top, gg2list(y))
  return(x)
}

##' @inheritParams atleft
##' @rdname loc
##' @export
atbottom <- function(x, y) {
  x@bottom <- c(x@bottom, gg2list(y))
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


##' Transfer to a list
##'
##' @title Transfer to list.
##' @param x A \code{ggplot}/\code{gtable} object.
##' @return A list.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom ggplot2 ggplotGrob
##' @keywords internal
##' 
gg2list <- function(x) {
  if (inherits(x, 'ggplot')) {
    x <- list(ggplotGrob(x))
  }
  else if (inherits(x, 'gtable')){
    x <- list(x)
  }
  else {}

  return(x)
}


##' Convert to the core of a \code{gmat} class object.
##'
##' @title Convert to the core
##' @param x A \code{ggplot} or \code{gtable} object.
##' @return A \code{pmat} object.
##' @examples
##' require('ggplot2')
##' require('magrittr')
##'
##' qplot(1,1) %>% ascore
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##' 
ascore <- function(x) {
  return(new('gmat', core = gg2list(x)))
}
