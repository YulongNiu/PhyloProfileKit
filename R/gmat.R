##' @include AllClasses.R AllGenerics.R operators.R
NULL

plot.gmat <- function(x, ...) {
  nleft <- ncol(x@left)
  nright <- ncol(x@right)
  ntop <- ncol(x@top)
  nbottom <- ncol(x@bottom)




}


AddEmpty <- function(x, conNum, rowNum, rev = FALSE) {

  ## plot empty block
  eObj <- gg2t(pp_empty(colour = 'white'))
}

