##' @include AllClasses.R
NULL


##' Show method for \code{PP} objects
##'
##' @title show methods
##' @param object A \code{PP} object.
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
            d <- object@.Data

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profile"\n')
            cat('class: ', class(object), '\n')
            if (is.integer(object)) {
              cat('type: "binning"', '\n')
            } else {
              cat('type: "continuous"', '\n')
            }
            cat('#species: ', ncol(d), '\n')
            cat('#gene/protein: ', nrow(d), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(d)

          })


##' The constructor the PP class
##'
##' Construct a PP object from a matrix.
##'
##' @title Constructor
##' @param value A numeric matrix.
##' @return A PP object.
##' @examples
##' require("magrittr")
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
##' @importFrom magrittr %>%
##' @export
##'
PP <- function(value) {

  if (class(value) == 'matrix') {

    colSize <- ncol(value)
    rowSize <- nrow(value)

    ## check 0 row or 0 columns
    if (colSize == 0 ||
        rowSize == 0) {
      value %>% as('PP') %>% return
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
    value %>% as('PP') %>% return
  } else {
    value %>% return
  }
}



##' Select and replace the data matrix of a PP object
##'
##' \code{PPData(x)}: Extract the data matrix from a PP object.
##'
##' \code{PPData(x) <- value}: Replace the data matrix of a PP object.
##'
##' @title Select and replace PP data
##' @inheritParams PPData
##' @return A numeric matrix.
##'
##' \code{PPData(x)}: A numeric matrix.
##'
##' \code{PPData(x) <- value}: An update PP object.
##'
##' @examples
##' require(magrittr)
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
##' @importFrom magrittr %>%
##' @rdname PPData-methods
##' @exportMethod PPData
##'
setMethod(f = 'PPData',
          signature = c(x = 'PP'),
          definition = function(x, ...) {
            x@.Data %>% return
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


setMethod(f = '[',
          signature = c(x = 'PP'),
          definition = function(x, i, j, ..., drop) {
            x@.Data <- callNextMethod()
            x %T>% validObject %>% return
          })


## select and replace

## cbind

## rbind

