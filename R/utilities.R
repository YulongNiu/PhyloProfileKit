##' Additional utilities
##'
##' \code{isBinMat_internal()}: Whether a matrix is binning.
##'
##' \code{valiMat_internal()}: Validate a numeric matrix with row names and column names.
##' 
##' @title utilities
##' @param x A matrix.
##' @return
##'
##' \code{isBinMat_internal()}: Logic.
##'
##' \code{valiMat_internal()}: \code{TRUE} or a warning message.
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


##' @param x A numeric matrix.
##' @param className A character string indicating the class names.
##' @rdname utilities
##' @keywords internal
##' 
valiMat_internal <- function(x, className) {

  warnNumMat <- 'The %class% (.Data slot) should be an integer matrix of a numeric matrix.'
  warnMatName <- 'The %class% (.Data slot) needs rownames or colnames.'

  ## 1. validate numeric matrix
  if (!(is.numeric(x) &&
        is.matrix(x))) {
    warn <- sub('%class%', className, warnNumMat)
    return(warn)
  } else {}

  ## 2. validate rownames and colnames
  if (ncol(x) == 0 ||
      nrow(x) == 0){
    return(TRUE)
  } else {
    if (is.null(rownames(x)) ||
        is.null(rownames(x))) {
      warn <- sub('%class%', className, warnMatName)
      return(warn)
    } else {
      return(TRUE)
    }
  }
}
