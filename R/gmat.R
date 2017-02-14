##' @include AllClasses.R AllGenerics.R operators.R
NULL

plot.gmat <- function(x, ...) {
  nleft <- length(x@left)
  nright <- length(x@right)
  ntop <- length(x@top)
  nbottom <- length(x@bottom)




}


AddEmpty <- function(x, leftN, rightN, coreN, reverse = FALSE) {

  if(reverse) {
    x <- rev(x)
  } else {}

  ## plot empty block
  eObj <- ggplotGrob(pp_empty(colour = 'white'))

  rowN <- length(x)
  colN <- leftN + rightN + coreN

  l <- lapply(1:(colN * rowN), function(x){eObj})
  insertIdx <- seq(1, colN * rowN, by = colN) + leftN
  l[insertIdx] <- x

  return(l)
}

