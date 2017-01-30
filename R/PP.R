##' @include AllClasses.R AllGenerics.R utilities.R
NULL


##' Show method for \code{PP} and \code{PPIdx} objects
##'
##' @title Show methods
##' @param object A \code{PP}/\code{PPIdx}/\code{PPResult}/\code{PPTreeIdx} object.
##' @return Show messages.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
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

            str(p)

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
##' @author Yulong Niu \email{niuylscu@@gmail.com}
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
##' @return A numeric matrix.
##'
##' \code{PPData(x)}: A numeric matrix.
##'
##' \code{PPData(x) <- value}: An update PP object.
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
##' rmat <- matrix(sample(0:1, 20 * 30, replace = TRUE),
##'                ncol = 30,
##'                dimnames = list(paste0('pro', 1:20),
##'                                paste0('spe', 1:30)))
##' PPData(ppContinuous)  <- rmat
##' ## type is changed to "binning"
##' ppContinuous
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
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
          signature = c(x = 'PP'),
          definition = function(x, ..., value) {
            x@.Data <- value
            x %T>% validObject %>% return
          })



##' Select or replace parts of a \code{PP} object
##'
##' \code{x[i, j, ..., drop]}: Select parts of a \code{PP} object.
##'
##' \code{x[i, j, ..., drop] <- value}: Replace parts of a \code{PP} object.
##'
##' @title Select or replace PP objects
##' @param x A \code{PP} object.
##' @param i,j,... Indices.
##' @param drop Whether the result is coerced to the lowest possible dimension. Desalt set is \code{FALSE}.
##' @return A \code{PP} object.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
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
          signature = c(x = 'PP'),
          definition = function(x, i, j, ..., value) {
            PPData(x) <- callNextMethod()
            return(x)
          })



## TODO:
## cbind
## rbind
