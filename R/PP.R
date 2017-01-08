##' @include AllClasses.R show-methods.R
NULL

##' The constructor and accessors of the PP class
##'
##' PP(): Construct a PP object from a matrix.
##'
##' PPData(): Extract the data matrix from a PP object. See the following examples for replacing the data matrix.
##'
##' @title Constructor and acccessors
##' @param value A numeric matrix.
##' @return
##'
##' PP(): A PP object.
##'
##' PPData(): A numeric matrix.
##'
##' @examples
##'
##' ## construct a PP object without dimnames
##' ppBinning <- PP(matrix(sample(0:1, 10 * 20, replace = TRUE), ncol = 20))
##'
##' ## with dimnames
##' ppContinuous <- new('PP',
##'                     matrix(rnorm(10 * 20),
##'                            ncol = 20,
##'                            dimnames = list(paste0('protein', 1:10),
##'                                            paste0('spe', 1:20))))
##'
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
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname PP
##' @export
##'
PP <- function(value) {

  if (class(value) == 'matrix') {

    colSize <- ncol(value)
    rowSize <- nrow(value)

    ## check 0 row or 0 columns
    if (colSize == 0 ||
        rowSize == 0) {
      return(as(value, 'PP'))
    } else {}

    ## check colnames
    if (is.null(colnames(value))) {
      colnames(value) <- paste0('spe', seq_len(colSize))
    } else {}

    ## check rownames
    if (is.null(rownames(value))) {
      rownames(value) <- paste0('gene', seq_len(rowSize))
    } else{}

    ## transfer
    return(as(value, 'PP'))
  } else {
    return(value)
  }
}

##' @param ppObj A PP object.
##' @rdname PP
##' @export
##'
PPData <- function(ppObj) {
  if (is(ppObj, 'PP')) {
    return(ppObj@.Data)
  } else {
    return(ppObj)
  }
}


##' @inheritParams PP
##' @inheritParams PPData
##' @importFrom methods validObject
##' @rdname PP
##' @export
##'
'PPData<-' <- function(ppObj, value) {
  if (is(ppObj, 'PP')) {
    ppObj@.Data <- value
    validObject(ppObj)
    return(ppObj)
  } else {
    return(value)
  }
}



## select and replace

## cbind

## rbind



