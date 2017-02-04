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


##' @param x Object with at least profile and linkage indices.
##' @param n The number of CPUs or processors.
##' @param FUN Functions to process paired linkages. The \code{FUN} has two parameters, the first parameter is a list \code{list(f = f, t = t, uniID = i)} in which \code{f} and \code{t} is two numeric vector indicating two single profile and \code{uniID} is the current index for parallel programming, and the second parameter is \code{...} used for additional arguments.
##' @param ... Additional parameters passed to \code{FUN}.
##' @rdname Batch-methods
##' @export
##' 
setGeneric(name = 'Batch',
           def = function(x, FUN, ..., n){standardGeneric('Batch')})


##' @param x A \code{PPIdx}/\code{PPTreeIdx} object.
##' @param method A character string.
##' \itemize{
##'   \item \code{"SimCor"}: Person's correlation coefficient.
##'   \item \code{"SimJaccard"}: Jaccard similarity.
##'   \item \code{"SimMI"}: Mutual information.
##'   \item \code{"DistHamming"}: Hamming distance.
##'   \item \code{"DistEuclidean"}: Euclidean distance.
##' }
##'
##' \code{"SimCor"} and \code{"DistHamming"} are not recommended for binning profiles.
##' @param ... Additional parameters passed to \code{method}.
##' \itemize{
##'   \item \code{bin}: Integer. The number of breaks in the \code{"MI"} method for continuous profiles.
##' }
##' @inheritParams Batch
##' @rdname SimDist-methods
##' @export
##' 
setGeneric(name = 'SimDist',
           def = function(x, method, ..., n){standardGeneric('SimDist')})


##' @param x A \code{PP} object.
##' @inheritParams SimDist
##' @rdname ChooseSimDistFun-methods
##' @keywords internal
##' 
setGeneric(name = 'ChooseSimDistFun',
           def = function(x, method, ...){standardGeneric('ChooseSimDistFun')})



##' @param x \code{PPTreeIdx} object.
##' @param ... Additional parameters.
##' @inheritParams Batch
##' @rdname Dollo-methods
##' @export
##' 
setGeneric(name = 'Dollo',
           def = function(x, ..., n){standardGeneric('Dollo')})


##' @param ... Additional parameters
##' \itemize{
##'   \item \code{BayesTraitsPath}: Full path of the program \href{http://www.evolution.reading.ac.uk/BayesTraits.html}{BayesTraits}.
##'   \item \code{method}: Either "MCMC" or "ML" (Maximum Likelihood) method. "ML" is set as default.
##' }
##' Parameters if \code{method} is set as \code{"MCMC"}:
##' \itemize{
##'   \item \code{priorAll}: Prior distribution.
##'   \item \code{iterNum}:  Iteration number with default value 1010000.
##' }
##' @inheritParams Dollo
##' @rdname BayesTraits-methods
##' @export
##' 
setGeneric(name = 'BayesTraits',
           def = function(x, ..., n) {standardGeneric('BayesTraits')})



