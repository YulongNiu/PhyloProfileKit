##' @include AllClasses.R AllGenerics.R
NULL

##' Parallel framework for \code{PPIdx} and \code{PPTreeIdx}.
##'
##' \code{Batch()}: Parallel framework to analysis paired profiles in batch mode.
##'
##' @title Parallel framework
##' @inheritParams Batch
##' @return A numeric vector.
##' @examples
##' require('magrittr')
##'
##' ## Person correlation coefficient
##' ppBinIdx <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PP %>% PPIdx(1:3, 1:3)
##' testfun <- function(eachArg, ...) {sum(eachArg$f * eachArg$t)}
##' Batch(ppBinIdx, testfun, n = 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @rdname Batch-methods
##' @seealso Batch process of similarity/distance \code{\link{SimDist}}.
##' @exportMethod Batch
##'
setMethod(f = 'Batch',
          signature = c(x = 'PPIdx'),
          definition = function(x, FUN, ..., n = 1) {

            ## register multiple core
            registerDoParallel(cores = n)

            p <- PPData(x)
            ft <- x@idx
            ppiNames <- rownames(ft)
            ppiNum <- nrow(ft)

            batchVec <- foreach(i = 1:ppiNum, .combine = c) %dopar% {
              ## print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
              f <- p[ft[i, 1], ]
              t <- p[ft[i, 2], ]
              eachSD <- FUN(eachArg = list(f = f, t = t, uniID = i), ...)

              return(eachSD)
            }

            ## stop multiple core
            stopImplicitCluster()

            return(batchVec)
          })
