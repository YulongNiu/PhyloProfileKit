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
##' Batch(ppBinIdx, testfun)
##' Batch(ppBinIdx, testfun, n = 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname Batch-methods
##' @seealso Batch process of similarity/distance \code{\link{SimDist}}.
##' @exportMethod Batch
##'
setMethod(f = 'Batch',
          signature = c(x = 'PPIdx'),
          definition = function(x, FUN, ..., n = 1) {

            p <- PPData(x)
            idx <- IdxData(x)

            batchVec <- BatchCore(p = p,
                                  idx = idx,
                                  FUN = FUN,
                                  ...,
                                  n = n)
            return(batchVec)
          })


##' Parallel core for \code{numeric matrix} and \code{big.matrix}.
##'
##' Parallel core for the index data.
##'
##' @title Parallel core
##' @inheritParams BatchCore
##' @return A \code{numeric matrix} object.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @rdname BatchCore-methods
##' @keywords internal
##'
setMethod(f = 'BatchCore',
          signature = c(idx = 'matrix'),
          definition = function(p, idx, FUN, ..., n = 1) {

            ## register multiple core
            registerDoParallel(cores = n)

            itx <- iter(idx, by = 'row')

            batchVec <- foreach(i = itx, .combine = c) %dopar% {
              f <- p[i[1], ]
              t <- p[i[2], ]
              eachVal <- FUN(eachArg = list(f = f, t = t, uniID = i), ...)
              return(eachVal)
            }

            ## stop multiple cores
            stopImplicitCluster()

            print('It is using numeric matrix.\n')

            return(batchVec)
          })


##' @inheritParams BatchCore
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators icount
##' @importFrom magrittr %>%
##' @rdname BatchCore-methods
##' @keywords internal
setMethod(f = 'BatchCore',
          signature = c(idx = 'big.matrix'),
          definition = function(p, idx, FUN, ..., n = 1) {

            ## register multiple core
            registerDoParallel(cores = n)

            itx <- idx %>% nrow %>% icount

            batchVec <- foreach(i = itx, .combine = c) %dopar% {
              f <- p[idx[i, 1], ]
              t <- p[idx[i, 2], ]
              eachVal <- FUN(eachArg = list(f = f, t = t, uniID = i), ...)
              return(eachVal)
            }

            ## stop multiple cores
            stopImplicitCluster()

            print('It is using big.matrix.\n')

            return(batchVec)

          })
