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


##' @param p A \code{PP} object.
##' @param x A matrix with two columns or a numeric vector. Proteins not in \code{p} (profile) are removed. The numeric vector indicates the indices of proteins.
##' \itemize{
##'   \item a \code{character matrix}: Each row containing the protein names is a linkage.
##'   \item a \code{numeric matrix} or \code{big.matrix}: Each row containing the protein indices is a linkage.
##'   \item a \code{numeric vector}: Each element is the index of a interested protein.
##' }
##' @param ... Additional parameters if \code{x} is a numeric vector.
##' \itemize{
##'   \item \code{y}: Another numeric vector used to generate paired linkages with \code{x}. Each element of \code{x} should be in \code{y}.
##'   \item \code{self}: Whether include self pairs, and default set is \code{FALSE}.
##'   \item \code{bidirect}: Whether to include two directions, and default set is \code{FALSE}.
##' }
##' @rdname Idx-methods
##' @export
##'
setGeneric(name = 'Idx',
           def = function(p, x, ...){standardGeneric('Idx')})


##' @param x A \code{PPIdx}/\code{PPTreeIdx} object.
##' @param ... Additional parameters.
##' @rdname IdxData-methods
##' @export
##'
setGeneric(name = 'IdxData',
           def = function(x, ...){standardGeneric('IdxData')})


##' @param value A numeric matrix.
##' @inheritParams IdxData
##' @rdname IdxData-methods
##' @export
##'
setGeneric(name = 'IdxData<-',
           def = function(x, ..., value){standardGeneric('IdxData<-')})


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


##' @param x A \code{PPIdx}/\code{PPTreeIdx} object.
##' @param n The number of CPUs or processors.
##' @param method Functions to process paired linkages.
##' \itemize{
##'   \item \code{SimCor}: Person's correlation coefficient.
##'   \item \code{SimJaccard}: Jaccard similarity.
##'   \item \code{SimMI}: Mutual information.
##'   \item \code{DistHamming}: Hamming distance.
##'   \item \code{DistManhattan}: Manhattan distance.
##'   \item \code{DistEuclidean}: Euclidean distance.
##'   \item \code{DistMinkowski}: Minkowski distance.
##'   \item \code{custom}: Custom similarity/distance functions.
##' }
##' The method \code{DistHamming} and \code{DistManhattan} generated the same distance value for binary phylogenetic profile. The \code{DistManhattan} is recommended for the continuous phylgenetic profile.
##' @param ... Additional parameters passed to \code{method}.
##' \itemize{
##'   \item \code{bin}: The bin size for the \code{SimMiConti} method, and the default value is 10.
##'    \item \code{p}: A \code{integer} indicating the \code{p} parameter for the \code{DistMinkowski} method, and the default value is 3.
##'    \item \code{func}: The pointer to the custom function.
##' }
##' @rdname Batch-methods
##' @export
##'
setGeneric(name = 'Batch',
           def = function(x, method, ..., n = 1){standardGeneric('Batch')})


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
           def = function(x, method, ..., n = 1){standardGeneric('SimDist')})


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
           def = function(x, ..., n = 1){standardGeneric('Dollo')})


##' @param ... Additional parameters
##' \itemize{
##'   \item \code{BayesTraitsPath}: Full path of the program \href{http://www.evolution.reading.ac.uk/BayesTraits.html}{BayesTraits}.
##'   \item \code{method}: Either "MCMC" or "ML" (Maximum Likelihood) method. "ML" is set as default.
##' }
##' Parameters only used in \code{MCMC}:
##' \itemize{
##'   \item \code{priorAll}: Prior distribution.
##'   \item \code{iterNum}:  Iteration number with default value 1010000.
##' }
##' @inheritParams Dollo
##' @rdname BayesTraits-methods
##' @export
##' 
setGeneric(name = 'BayesTraits',
           def = function(x, ..., n = 1) {standardGeneric('BayesTraits')})


##' @param x \code{PP}/\code{PPTreeIdx} object.
##' @param method Distance method used in the \code{dist()}.
##' @param ... Additional parameters.
##' \itemize{
##'   \item \code{proSize}: A numeric value. Size of protein names.
##'   \item \code{proGroup}: A factor indicating protein groups.
##'   \item \code{proGroupCol}: A named character vector indicating protein group colour.
##'   \item \code{speGroup}: A factor indicating species groups.
##'   \item \code{speGroup}: A named character vector indicating species group colour.
##' }
##' @rdname plotprofile-methods
##' @export
##' 
setGeneric(name = 'plotprofile',
           def = function(x, method, ...) {standardGeneric('plotprofile')})
