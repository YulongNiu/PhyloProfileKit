##' Additional utilities
##'
##' \code{isBinMat_internal()}: Whether a matrix is binning.
##'
##' \code{valiMatNames_internal()}: Validate a numeric matrix with row names and column names.
##' 
##' @title utilities
##' @param x A matrix.
##' @return
##'
##' \code{isBinMat_internal()}: Logic.
##'
##' \code{valiMatNames_internal()}: \code{TRUE} or a warning message.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
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
    return(identical(uniqVec, 0:1))
  }
  else {
    return(FALSE)
  }
}


##' @param x A numeric matrix.
##' @param warnName A character string indicating the class names.
##' @importFrom bigmemory is.big.matrix
##' @rdname utilities
##' @keywords internal
##' 
valiMatNames_internal <- function(x, warnName) {

  warnNumMat <- '%name% should be a matrix or big.matrix.'
  warnMatName <- '%name% needs rownames or colnames.'

  ## 1. validate numeric matrix
  if (!(is.matrix(x) ||
        is.big.matrix(x))) {
    warn <- sub('%name%', warnName, warnNumMat)
    return(warn)
  } else {}
  ## 2. validate rownames and colnames
  if (ncol(x) == 0 ||
      nrow(x) == 0){
    return(TRUE)
  } else {
    if (is.null(rownames(x)) ||
        is.null(rownames(x))) {
      warn <- sub('%name%', warnName, warnMatName)
      return(warn)
    } else {
      return(TRUE)
    }
  }
}

##' @param m A numeric matrix.
##' @param t A phylo tree.
##' @rdname utilities
##' @keywords internal
##' 
validTreeMat_internal <- function(m, t) {
  if (!all(colnames(m) == t$tip.label)) {
    return('rownames of profile should be in same order with tree tips.')
  } else {
    return(TRUE)
  }
}

##' Combination functions
##'
##' \code{comb2_interal()}: Pair combination of the input vector \code{x}.
##'
##' \code{combWhole_internal()}: Pair combination of the input vector \code{x} and \code{y}.
##'
##' \code{combSelf_internal()}: Self pairs.
##'
##' @title Combination
##' @param x Vector.
##' @return
##'
##' \code{comb2_internal()}/\code{combWhole_internal()}/\code{combSelf_internal()}: Matrix with 2 columns
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname comb
##' @keywords internal
##' 
comb2_internal <- function(x) {

  l <- length(x)

  if (l < 2) {
    m <- matrix(0, nrow = 0, ncol = 2)
  } else {
    fV <- rep(x[1:(l - 1)], ((l - 1):1))
    tVIdx <- unlist(mapply(seq, 2:l, l))
    tV <- x[tVIdx]

    m <- cbind(fV, tV)
  }

  return(m)
}


##' @param y Another vector, and every element of \code{x} should be in \code{y}.
##' @param self Whether include self pairs.
##' @param bidirect Whether to include two directions.
##' @rdname comb
##' @keywords internal
##' 
combWhole_internal <- function(x, y, self = FALSE, bidirect = FALSE) {

  ## check x
  ## first unique x
  x <- unique(x)
  x <- x[x %in% y]

  l <- length(x)

  if (l >= length(y)) {
    m <- comb2_internal(x)
  } else {
    yLeft <- y[!(y %in% x)]
    m <- cbind(rep(x, each = length(yLeft)),
               rep(yLeft, l))
    m <- rbind(m,
               comb2_internal(x))
  }

  if (bidirect) {
    m <- rbind(m, m[, 2:1])
  } else {}

   if (self) {
    m <- rbind(m, combSelf_internal(x))
   } else {}

  colnames(m) <- NULL

  return(m)
}


##' @inheritParams comb2_internal()
##' @rdname comb
##' @keywords internal
combSelf_internal <- function(x) {
  m <- cbind(x, x)
  return(m)
}
