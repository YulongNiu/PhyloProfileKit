##' @param x A PP object.
##' @param ... Additional parameters
##' @rdname PPData-methods
##' @export
##' 
setGeneric(name = 'PPData',
           def = function(x, ...){standardGeneric('PPData')})


##' @param value A numeric matrix.
##' @inheritParams PPData
##' @rdname PPData-methods
##' @export
##'
setGeneric(name = 'PPData<-',
           def = function(x, ..., value){standardGeneric('PPData<-')})
