##' @include AllClasses.R AllGenerics.R utilities.R
NULL


##' Show method for \code{PP} and \code{PPIdx} objects
##'
##' @title Show methods
##' @param object A \code{PP}/\code{PPIdx}/\code{PPTree}/\code{PPTreeIdx}/code{PPResult} object.
##' @return Show messages.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom utils str
##' @importFrom methods show
##' @rdname show-methods
##' @exportMethod show
##'
setMethod(f = 'show',
          signature = 'PP',
          definition = function(object){

            p <- object@.Data

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profile"\n')
            cat('class: ', class(object), '\n')
            if (isBinMat_internal(p)) {
              cat('profile: binning', '\n')
            } else {
              cat('profile: continuous', '\n')
            }
            cat('#species: ', ncol(p), '\n')
            cat('#proteins: ', nrow(p), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(object)

          })


##' The constructor the \code{PP} class
##'
##' Construct a \code{PP} object from a numeric matrix.
##'
##' @title Constructor of \code{PP}
##' @param x A numeric matrix.
##' @return A \code{PP} object.
##' @examples
##' require('magrittr')
##'
##' ## construct a PP object without dimnames
##' ppBinning <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PP
##'
##' ## with dimnames
##' ppContinuous <- matrix(rnorm(10 * 20),
##'                        ncol = 20,
##'                        dimnames = list(paste0('protein', 1:10),
##'                                        paste0('spe', 1:20))) %>% PP
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom methods as
##' @export
##'
PP <- function(x) {

  if (is.matrix(x) &&
      is.numeric(x)) {

    colSize <- ncol(x)
    rowSize <- nrow(x)

    ## check 0 row or 0 columns
    if (colSize == 0 ||
        rowSize == 0) {
      return(as(x, 'PP'))
    } else {}

    ## check colnames
    if (is.null(colnames(x))) {
      colnames(x) <- paste0('spe', seq_len(colSize))
    } else {}

    ## check rownames
    if (is.null(rownames(x))) {
      rownames(x) <- paste0('protein', seq_len(rowSize))
    } else{}

    ## transfer
    return(as(x, 'PP'))
  } else {
    return(x)
  }
}



##' Select and replace the data matrix of a \code{PP} object
##'
##' \code{PPData(x)}: Extract the data matrix from a \code{PP} object.
##'
##' \code{PPData(x) <- value}: Replace the data matrix of a \code{PP} object.
##'
##' @title Select and replace \code{PP} data
##' @inheritParams PPData
##' @return
##'
##' \code{PPData(x)}: A numeric matrix.
##'
##' \code{PPData(x) <- value}: An update \code{PP} object.
##'
##' @examples
##' require('magrittr')
##' ppContinuous <- matrix(rnorm(10 * 20),
##'                        ncol = 20,
##'                        dimnames = list(paste0('protein', 1:10),
##'                                        paste0('spe', 1:20))) %>% PP
##'
##' ## extract data matrix
##' PPData(ppContinuous)
##'
##' ## replace whole data matrix
##' pmat <- matrix(sample(0:1, 20 * 30, replace = TRUE),
##'                ncol = 30,
##'                dimnames = list(paste0('pro', 1:20),
##'                                paste0('spe', 1:30)))
##' PPData(ppContinuous)  <- pmat
##' ## type is changed to "binning"
##' ppContinuous
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname PPData-methods
##' @exportMethod PPData
##'
setMethod(f = 'PPData',
          signature = c(x = 'PP'),
          definition = function(x, ...) {
            return(x@.Data)
          })


##' @inheritParams PPData
##' @importFrom methods validObject
##' @importFrom magrittr %>% %T>%
##' @rdname PPData-methods
##' @exportMethod PPData<-
##'
setMethod(f = 'PPData<-',
          signature = c(x = 'PP', value = 'matrix'),
          definition = function(x, ..., value) {
            x@.Data <- value
            x %T>% validObject %>% return
          })



##' Select or replace parts of a \code{PP}/\code{PPResult} object
##'
##' \code{x[i, j, ..., drop]}: Select parts of a \code{PP}/\code{PPResult} object.
##'
##' \code{x[i, j, ..., drop] <- value}: Replace parts of a \code{PP} object.
##'
##' @title Select or replace \code{PP}/\code{PPResult} objects
##' @param x A \code{PP}/\code{PPResult} object.
##' @param i,j,... Indices.
##' @param drop Whether the result is coerced to the lowest possible dimension. Desalt set is \code{FALSE}.
##' @return A \code{PP}/\code{PPResult} object.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom methods callNextMethod
##' @rdname select-methods
##' @exportMethod [
##'
setMethod(f = '[',
          signature = c(x = 'PP'),
          definition = function(x, i, j, ..., drop = FALSE) {
            PPData(x) <- callNextMethod()
            return(x)
          })


##' @param value A number, numeric vector or numeric matrix.
##' @importFrom methods callNextMethod
##' @rdname select-methods
##' @exportMethod [<-
##'
setMethod(f = '[<-',
          signature = c(x = 'PP', value = 'numeric'),
          definition = function(x, i, j, ..., value) {
            PPData(x) <- callNextMethod()
            return(x)
          })


##' Row bind or column bind of a \code{PP} object
##'
##' \code{rbind(..., deparse.level = 1)}: Row bind a \code{PP} object with another \code{PP} object or a named numeric matrix.
##'
##' \code{cbind(..., deparse.level = 1)}: Column bind a \code{PP} object with another \code{PP} object or a named numeric matrix.
##'
##' @title rbind or cbind PP objects
##' @param x A \code{PP} object .
##' @param y A \code{PP} object or a named numeric matrix.
##' @param ... Additional parameters.
##' @return A \code{PP} object.
##' @examples
##' require('magrittr')
##' ppC <- matrix(rnorm(10 * 20),
##'               ncol = 20,
##'               dimnames = list(paste0('protein', 1:10),
##'                               paste0('spe', 1:20))) %>% PP
##'
##' ppM1 <- matrix(rnorm(10 * 20),
##'               ncol = 20,
##'               dimnames = list(paste0('protein', 1:10),
##'                               paste0('spe', 21:40)))
##'
##' ppM2 <- matrix(rnorm(10 * 20),
##'               ncol = 20,
##'               dimnames = list(paste0('protein', 11:20),
##'                               paste0('spe', 1:20)))
##'
##' cbind(ppC, ppM1)
##' cbind(ppC, PP(ppM1))
##' rbind(ppC, ppM2)
##' rbind(ppC, PP(ppM2))
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom methods rbind2 is
##' @rdname bind-methods
##' @exportMethod rbind2
##'
setMethod(f = 'rbind2',
          signature = c(x = 'PP'),
          definition = function(x, y, ...) {
            if (is(y, 'PP')) {
              PPData(x) <- rbind2(PPData(x), PPData(y))
            }
            else if (is(y, 'matrix')) {
              PPData(x) <- rbind2(PPData(x), y)
            }
            else {}

            return(x)
          })

##' @importFrom methods cbind2 is
##' @rdname bind-methods
##' @exportMethod cbind2
##'
setMethod(f = 'cbind2',
          signature = c(x = 'PP'),
          definition = function(x, y, ...) {
            if (is(y, 'PP')) {
              PPData(x) <- cbind2(PPData(x), PPData(y))
            }
            else if (is(y, 'matrix')) {
              PPData(x) <- cbind2(PPData(x), y)
            }
            else {}

            return(x)
          })

