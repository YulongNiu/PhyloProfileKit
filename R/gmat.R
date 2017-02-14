##' @include AllClasses.R AllGenerics.R

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

            leftN <- length(x@left)
            rightN <- length(x@right)
            topN <- length(x@top)
            bottomN <- length(x@bottom)
            coreN <- length(x@core)

            rowN <- topN + coreN + bottomN
            colN <- leftN + coreN + rightN

            pList <- c(AddEmpty(x@top,
                                leftN,
                                rightN,
                                coreN,
                                reverse = TRUE),
                       c(rev(x@left),
                         x@core,
                         x@right),
                       AddEmpty(x@bottom,
                                leftN,
                                rightN,
                                coreN,
                                reverse = TRUE))

            pObj <- grid.arrange(grobs = pList, ncol = colN, ...)

            return(pObj)
          })

##' Combine plot elements with empty blocks
##'
##' @title Combine empty blocks
##' @param x A list.
##' @param leftN left number.
##' @param rightN right number.
##' @param coreN center number.
##' @param reverse Reverse \code{x} or not.
##' @return A list.
##' @importFrom ggplot2 ggplotGrob
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##' 
AddEmpty <- function(x, leftN, rightN, coreN, reverse = FALSE) {

  if (length(x) == 0) {
    return(x)
  } else {}

  if (reverse) {
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

## library(PhyloProfile)
## library(ggplot2)
## library(gridExtra)
## library(magrittr)

## tmp1 <- matrix(sample(0:1, 10 * 40, replace = TRUE), ncol = 40)
## rownames(tmp1) <- paste0('protein', 1:10)
## colnames(tmp1) <- paste0('spe', 1:40)

## tmp2 <- pp_profile(tmp1) %>% ascore %@<% (pp_text(rownames(tmp1)) %@+% pp_tile(rownames(tmp1))) %@^% pp_tile(colnames(tmp1)) %@v% pp_text(colnames(tmp1))
