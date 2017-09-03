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
            }
            else if (is.numeric(x)) {
              x <- x
            } else {}

            return(x)
          })


##' @inheritParams Idx
##' @importFrom magrittr %<>%
##' @rdname Idx-methods
##' @exportMethod Idx
##'
setMethod(f = 'Idx',
          signature = c(p = 'PP', x = 'big.matrix'),
          definition = function(p, x, ...) {
            return(x)
          })


##' @inheritParams Idx
##' @importFrom magrittr %<>%
##' @rdname Idx-methods
##' @exportMethod Idx
##'
setMethod(f = 'Idx',
          signature = c(p = 'PP', x = 'vector'),
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
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom methods new
##' @importFrom bigmemory as.big.matrix
##' @seealso PPTreeIdx
##' @export
##'
PPIdx <- function(p, x, ...) {

  x <- Idx(p, x, ...)

  colSize <- ncol(x)
  rowSize <- nrow(x)

  ## check 0 row or 0 columns
  if (colSize == 0 ||
      rowSize == 0) {
    return(new('PPIdx', p, idx = x))
  } else {}

  return(new('PPIdx', p, idx = x))
}


##' Select and replace the indices of a \code{PPIdx}/\code{PPTreeIdx} object
##'
##' \code{IdxData(x)}: Extract the indices from a \code{PPIdx}/\code{PPTreeIdx} object.
##'
##' \code{IdxData(x) <- value}: Replace the indices of a \code{PPIdx}/\code{PPTreeIdx} object.
##'
##' @title Select and replace \code{PPIdx}/\code{PPTreeIdx} data
##' @inheritParams IdxData
##' @return
##'
##' \code{PPData(x)}: A numeric matrix.
##'
##' \code{PPData(x) <- value}: An update \code{PPIdx}/\code{PPTreeIdx} object.
##'
##' @examples
##' require('magrittr')
##' ppContinuous <- matrix(rnorm(10 * 20),
##'                        ncol = 20,
##'                        dimnames = list(paste0('protein', 1:10),
##'                                        paste0('spe', 1:20))) %>% PP
##' ppi <- PPIdx(ppContinuous, 1:3, 1:3)
##'
##' ## extract indices
##' IdxData(ppi)
##'
##' ## replace indices
##' imat <- matrix(sample(1:10, 2 * 3, replace = FALSE),
##'                ncol = 2,
##'                dimnames = list(paste0('link', 1:3),
##'                                c('from', 'to')))
##' IdxData(ppi)  <- imat
##' ## type is changed to "binning"
##' ppContinuous
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname IdxData-methods
##' @exportMethod IdxData
##'
setMethod(f = 'IdxData',
          signature = c(x = 'PPIdx'),
          definition = function(x, ...) {
            return(x@idx)
          })


##' @inheritParams IdxData
##' @importFrom methods validObject
##' @importFrom magrittr %>% %T>%
##' @rdname IdxData-methods
##' @exportMethod IdxData<-
##'
setMethod(f = 'IdxData<-',
          signature = c(x = 'PPIdx', value = 'matrix'),
          definition = function(x, ..., value) {
            x@idx <- value
            x %T>% validObject %>% return
          })

