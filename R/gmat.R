##' @include AllClasses.R AllGenerics.R operators.R
NULL

plot.gmat <- function(x, ...) {
  nleft <- length(x@left)
  nright <- length(x@right)
  ntop <- length(x@top)
  nbottom <- length(x@bottom)




}


AddEmpty <- function(x, leftN, rightN, reverse = FALSE) {

  if(reverse = TRUE) {
    tmp1 <- rev(tmp1)
  } else {}

  ## plot empty block
  eObj <- ggplotGrob(pp_empty(colour = 'white'))

  l <- lapply(1:(rowN, colN), function(x){eObj})


}

