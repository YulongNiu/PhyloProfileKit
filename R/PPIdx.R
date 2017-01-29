##' @include PP.R
NULL

##' @inheritParams show
##' @importFrom utils str
##' @importFrom methods show
##' @rdname show-methods
##' @exportMethod show
##
setMethod(f = 'show',
          signature = 'PPIdx',
          definition = function(object){

            p <- object@.Data
            idx <- object@idx

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profile with linkage indices"\n')
            cat('class: ', class(object), '\n')
            if (isBinMat_internal(p)) {
              cat('profile: "binning"', '\n')
            } else {
              cat('profile: "continuous"', '\n')
            }
            cat('#species: ', ncol(p), '\n')
            cat('#proteins: ', nrow(p), '\n')
            cat('#linkages: ', nrow(idx), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(p)
            str(idx)

            })


##' The constructor the \code{PPIdx} class
##'
##' Construct a \code{PP} object.
##'
##' @title Constructor of \code{PPIdx}
##' @param p A \code{PP} object
##' @param x A character matrix with two columns or a numeric vector. Proteins not in \code{p} (profile) are removed. The numeric vector indicates the indices of proteins.
##' @param bigmat Whether store the indices as a big matrix. Set it as \code{TRUE} if the number of index is large.
##' @param ... Additional parameters if \code{x} is a numeric vector.
##' \itemize{
##'   \item \code{y}: Another numeric vector used to generate paired linkages with \code{x}. Every element of \code{x} should be in \code{y}.
##'   \item \code{self}: Whether include self pairs.
##'   \item \code{bidirect}: Whether to include two directions.
##' }
##' @return A \code{PPIdx} object.
##' @examples
##' require('magrittr')
##'
##' ppBinning <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PP
##'
##' ## pre-built linkages
##' linkM <- sample(1:30, 20 * 2, replace = TRUE) %>% paste0('protein', .) %>% matrix(ncol = 2)
##' PPIdx(ppBinning, linkM)
##'
##' ## within top 3 proteins
##' PPIdx(ppBinning, 1:3, 1:3)
##' ## top 3 proteins with whole profiles
##' PPIdx(ppBinning, 1:3, 1:nrow(ppBinning))
##' ## with self linkages
##' PPIdx(ppBinning, 1:3, 1:3, self = TRUE)
##' ## with bidirectional linkages
##' PPIdx(ppBinning, 1:3, 1:3, self = TRUE, bidirect = TRUE)
##'
##' ## bigmatrix
##' PPIdx(ppBinning, 1:10, 1:nrow(ppBinning), bigmat = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom methods new
##' @importFrom bigmemory as.big.matrix
##' @export
##' 
PPIdx <- function(p, x, ..., bigmat = FALSE) {

  ## whole proteins
  wp <- rownames(p)

  if (is.matrix(x) &&
      is.character(x)) {
    ## x is a character matrix
    ## check linkage proteins in wp, NA if not in wp
    x %<>% c %>% match(wp) %>% matrix(ncol = 2, dimnames = dimnames(x))
    hasLogic <- !(is.na(x[, 1]) | is.na(x[, 2]))
    x <- x[hasLogic, , drop = FALSE]
  }
  else if (is.numeric(x) &&
           is.null(dim(x))) {
    ## x is a character vector
    x %<>% combWhole_internal(...)
  }
  else {
    return(x)
  }

  colSize <- ncol(x)
  rowSize <- nrow(x)

  ## check 0 row or 0 columns
  if (colSize == 0 ||
      rowSize == 0) {
    return(new('PPIdx', p, idx = x))
  } else {}

  ## check colnames
  if (is.null(colnames(x))) {
    colnames(x) <- c('From', 'To')
  } else {}

  ## check rownames
  if (is.null(rownames(x))) {
    rownames(x) <- paste0('link', seq_len(rowSize))
  } else{}

  ## if big matrix
  if (bigmat) {
    x <- as.big.matrix(x)
  } else {}

  return(new('PPIdx', p, idx = x))
}

