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

            p <- PPData(object)
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

            str(object)

            })


##' Index linkages methods
##'
##' Index linkages from protein IDs or indices in batch.
##'
##' @title Index linkages
##' @inheritParams Idx
##' @return A \code{numeric matrix} object.
##' @examples
##' require('magrittr')
##'
##' ppBinning <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PP
##'
##' ## pre-built linkages
##' linkM <- sample(1:30, 20 * 2, replace = TRUE) %>% paste0('protein', .) %>% matrix(ncol = 2)
##' Idx(ppBinning, linkM)
##'
##' ## top 3 proteins with whole profiles
##' Idx(ppBinning, 1:3, 1:nrow(ppBinning))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname Idx-methods
##' @importFrom magrittr %>% %<>%
##' @exportMethod Idx
##'
setMethod(f = 'Idx',
          signature = c(p = 'PP', x = 'matrix'),
          definition = function(p, x, ...) {

            ## whole proteins
            wp <- rownames(p)

            if (is.character(x)) {
              ## x is a character matrix
              ## check linkage proteins in wp, NA if not in wp
              x %<>% c %>% match(wp) %>% matrix(ncol = 2, dimnames = dimnames(x))
              hasLogic <- !(is.na(x[, 1]) | is.na(x[, 2]))
              x <- x[hasLogic, , drop = FALSE]
            } else {}

            return(x)
          })


##' @inheritParams Idx
##' @importFrom magrittr %<>%
##' @rdname Idx-methods
##' @exportMethod Idx
##'
setMethod(f = 'Idx',
          signature = c(p = 'PP', x = 'numeric'),
          definition = function(p, x, ...) {
            if(is.null(dim(x))) {
              ## x is a numeric vector
              x %<>% combWhole_internal(...)
            } else {}
            return(x)
          })

##' The constructor the \code{PPIdx} class
##'
##' Construct a \code{PPIdx} object.
##'
##' @title Constructor of \code{PPIdx}
##' @param bigmat Whether store the indices as a big matrix. Set it as \code{TRUE} if the number of index is large.
##' @inheritParams Idx
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
##' @seealso PPTreeIdx
##' @export
##' 
PPIdx <- function(p, x, ..., bigmat = FALSE) {

  x <- Idx(p, x, ...)

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

