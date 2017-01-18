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
##' @param warnName A character string indicating the class names.
##' @rdname utilities
##' @keywords internal
##' 
valiMat_internal <- function(x, warnName) {

  warnNumMat <- 'The %name% should be an matrix.'
  warnMatName <- 'The %name% needs rownames or colnames.'

  ## 1. validate numeric matrix
  if (!is.matrix(x)) {
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

##' Combination functions
##'
##' \code{comb2_interal()}: Pair combination of the input vector \code{x}.
##'
##' \code{combWhole_internal()}: Pair combination of the input vector \code{x} and \code{y}.
##'
##' \code{combSelf_internal()}: Self pairs.
##'
##' @title Combination
##' @param x Vector with length > 0.
##' @param self Whether include self pairs.
##' @param bidirect Whether to include two directions.
##' @return
##'
##' \code{comb2_internal()}/\code{combWhole_internal()}/\code{combSelf_internal()}: Matrix with 2 columns
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname comb
##' @keywords internal
##' 
comb2_internal <- function(x, self = FALSE, bidirect = FALSE) {

  l <- length(x)

  if (l < 2) {
    if (!self) {
      return(matrix(0, 0, 0))
    } else {
      return(combSelf_internal(x))
    }
  } else {}

  fV <- rep(x[1:(l - 1)], ((l - 1):1))
  tVIdx <- unlist(mapply(seq, 2:l, l))
  tV <- x[tVIdx]

  m <- cbind(fV, tV)

  if (bidirect) {
    m <- rbind(m, m[, 2:1])
  } else {}

  if (self) {
    m <- rbind(m, combSelf_internal(x))
  } else ()

  return(m)
}


##' @param y The whole vector, and every element of \code{x} should be in \code{y}.
##' @inheritParams comb2_internal()
##' @rdname comb
##' @keywords internal
##' 
combWhole_internal <- function(x, y, self = FALSE, bidirect = FALSE) {

  l <- length(x)

  if (l < 2) {
    if (!self) {
      return(matrix(0, 0, 0))
    } else {
      return(combSelf_internal(x))
    }
  } else {}

  if (l == length(y)) {
    m <- comb2_interal(x)
  } else {
    yLeft <- y[!(y%in%x)]
    m <- cbind(rep(x, each = length(yLeft)),
               rep(yLeft, l))
    m <- rbind(m,
               comb2_interal(x))
  }

  if (bidirect) {
    m <- rbind(m, m[, 2:1])
  } else {}

   if (self) {
    m <- rbind(m, combSelf_internal(x))
   } else {}

  return(m)
}


##' @inheritParams comb2_internal()
##' @rdname comb
##' @keywords internal
combSelf_internal <- function(x) {
  m <- cbind(x, x)
  return(m)
}
