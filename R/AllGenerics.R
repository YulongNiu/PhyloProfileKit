##' Select and replace generic methods
##'
##' \code{PPData(x)}: Extract the data matrix from a PP object.
##'
##' \code{PPData(x) <- value}: Replace the data matrix of a PP object.
##'
##' @title Select and replace PP data
##' @param x A PP object.
##' @param ... Additional parameters
##' @return
##'
##' \code{PPData(x)}: A numeric matrix.
##'
##' \code{PPData(x) <- value}: An update PP object.
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname PPData-generic-methods
##' @export
##'
setGeneric(name = 'PPData',
           def = function(x, ...){standardGeneric('PPData')})


##' @param value A numeric matrix.
##' @inheritParams PPData
##' @rdname PPData-generic-methods
##' @export
##'
setGeneric(name = 'PPData<-',
           def = function(x, ..., value){standardGeneric('PPData<-')})
