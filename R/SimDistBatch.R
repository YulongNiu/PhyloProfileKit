##' @include utilities.R Batch.R
NULL

##' Similarity or distance.
##'
##' Similarity and distance of paired profiles.
##'
##' @inheritParams SimDist
##' @title Batch process of similarity and distance.
##' @return A numeric vector
##' @examples
##' require('magrittr')
##' data(fatp)
##'
##' ## Mutual information
##' f1 <- fatp$atpPhylo %>% PP %>% PPIdx(1:6, 1:6)
##' SimDist(f1, n = 2, 'SimMI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname SimDist-methods
##' @exportMethod SimDist
##'
setMethod(f = 'SimDist',
          signature = c(x = 'PPIdx'),
          definition = function(x, method, ..., n = 1) {

            if (isBinMat_internal(x@.Data)) {
              MI <- SimMIBin
            } else {
              MI <- SimMIConti
            }

            m <- switch(method,
                        SimCor = SimCor,
                        SimJaccard = SimJaccard,
                        SimMI = MI,
                        DistHamming = DistHamming,
                        DistEulidean = DistEuclidean)

            bv <- Batch(x = x, n = n, FUN = m, ...)

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = method)
            return(bvRes)
          })

