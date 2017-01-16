##' @param x A \code{PP} object.
##' @param ... Additional parameters.
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



##' @param x A \code{PP} object.
##' @param method A character string, "NPP" or "SVD". Method used to normalize raw bit score phylogenetic profile.
##' @param ... Additional parameters passed to normalization methods. Parameters shared in \code{NPP} and \code{SVD}:
##' \itemize{
##'   \item \code{bitCutoff}: Cutoff of minimum bit score.
##'   \item \code{bitReset}: Reset bit scores lower than the \code{bitCutoff}.

##'   \item \code{minConserve}: Minimum number of homologous. The proteins with homologous less than this value are discarded.
##' }
##' Parameters only used in \code{SVD}:
##' \itemize{
##'   \item \code{trimming}: Top percentages of the unitary matrix, for example 0.3 means top 30\% columns are kept.
##' }
##' @rdname Norm-methods
##' @export
##' 
setGeneric(name = 'Norm',
           def = function(x, method, ...){standardGeneric('Norm')})
