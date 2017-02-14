##' @include AllClasses.R AllGenerics.R operators.R
NULL

plot.gmat <- function(x, ...) {

  leftN <- length(x@left)
  rightN <- length(x@right)
  topN <- length(x@top)
  bottomN <- length(x@bottom)
  coreN <- length(x@core)

  rowN <- topN + coreN + bottomN
  colN <- leftN + coreN + rightN

  pList <- c(AddEmpty(x@top, leftN, rightN, coreN, reverse = TRUE),
            c(rev(x@left), x@core, x@right),
            AddEmpty(x@bottom, leftN, rightN, coreN, reverse = TRUE))

  pObj <- grid.arrange(grobs = pList, ncol = colN, ...)

  return(pObj)
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

