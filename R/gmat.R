##' @include AllClasses.R AllGenerics.R
NULL

##' Plot method for \code{gmat} objects
##'
##' @title Plot methods
##' @param x A \code{gmat}
##' @param y Not set ("missing").
##' @param ... Parameters passed to \code{grid.arrange()} in the gridExtra package.
##' @return A \code{gtable} object.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom graphics plot
##' @importFrom gridExtra grid.arrange
##' @rdname plot-methods
##' @exportMethod plot
##'
setMethod(f = 'plot',
          signature = c(x = 'gmat', y = 'missing'),
          definition = function(x, y, ...) {
            obj <- c(rev(x@top),
                     rev(x@left),
                     x@core,
                     x@right,
                     x@bottom)

            pObj <- grid.arrange(grobs = obj,
                                 layout_matrix = Lay(x),
                                 ...)

            return(pObj)
          })

##' Construct the layout matrix of a \code{gmat} object.
##'
##' @title Construct the layout matrix
##' @param x A \code{gmat} object.
##' @return A numeric matrix
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname lay
##' @keywords internal
##' 
Lay <- function(x) {
  leftN <- length(x@left)
  rightN <- length(x@right)
  coreN <- length(x@core)
  topN <- length(x@top)

  m <- rbind(LaySide(x, side = 'top'),
             seq_len(leftN + coreN + rightN) + topN,
             LaySide(x, side = 'bottom'))

  return(m)
}

##' @inheritParams Lay
##' @rdname lay
##' @keywords internal
##' 
LaySide <- function(x, side = 'top') {
  leftN <- length(x@left)
  rightN <- length(x@right)
  coreN <- length(x@core)
  topN <- length(x@top)
  bottomN <- length(x@bottom)

  if (side == 'top') {
    xn <- topN
    xIdx <- seq_len(xn)
  }
  else if (side == 'bottom') {
    xn <- bottomN
    if (length(xn) > 0) {
      xIdx <- seq_len(xn) + (topN + leftN + rightN + coreN)
    } else {
      xIdx <- seq_len(xn)
    }
  }
  else {
    stop('side must be "top" or "bottom".\n')
  }

  sideM <- matrix(ncol = (leftN + rightN + coreN),
                  nrow = xn)
  sideM[, (leftN + 1)] <- xIdx

  return(sideM)
}

