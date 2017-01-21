##' @include AllClasses.R AllGenerics.R utilities.R paraFrame.R
NULL

##' Similarity or distance.
##'
##' SimDist(): Similarity and distance of paired profiles.
##'
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
          definition = function(x, n = 1, method, ...) {

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
            return(bv)

          })

