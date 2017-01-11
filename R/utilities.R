##' Additional utilities
##'
##' \code{isBinMat_internal()}: Whether a matrix is binning.
##'
##' @title utilities
##' @param x A matrix.
##' @return
##'
##' \code{isBinMat_internal()}: Logic.
##'
##' 
##' @rdname utilities
##' @importFrom magrittr %>%
##' @keywords internal
##'
isBinMat_internal <- function(x) {
  uniqVec <- x %>% c %>% unique %>% sort
  lenV <- length(uniqVec)

  if (lenV == 1) {
    return(uniqVec == 0 || uniqVec == 1)
  }
  else if (lenV == 2) {
    return(all.equal(uniqVec, 0:1))
  }
  else {
    return(FALSE)
  }

}
