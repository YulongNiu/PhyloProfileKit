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
##' @rdname Batch-methods
##' @seealso Batch process of similarity/distance \code{\link{SimDist}}.
##' @exportMethod Batch
##'
setMethod(f = 'Batch',
          signature = c(x = 'PPIdx'),
          definition = function(x, FUN, ..., n = 1) {


            p <- PPData(x)
            idx <- x@idx

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
##' @rdname BatchCore-methods
##' @keywords internal
##'
setMethod(f = 'BatchCore',
          signature = c(p = 'matrix'),
          definition = function(p, idx, FUN, ..., n = 1) {

            ## register multiple core
            registerDoParallel(cores = n)

            ppiNum <- nrow(idx)

            batchVec <- foreach(i = 1:ppiNum, .combine = c) %dopar% {
              ## print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
              f <- p[idx[i, 1], ]
              t <- p[idx[i, 2], ]
              eachVal <- FUN(eachArg = list(f = f, t = t, uniID = i), ...)

              gc()

              return(eachVal)
            }

            ## stop multiple cores
            stopImplicitCluster()

            return(batchVec)
          })


##' @inheritParams BatchCore
##' @importFrom bigmemory describe attach.big.matrix
##' @importFrom parallel makeCluster clusterCall parSapply stopCluster
##' @importFrom doParallel registerDoParallel
##' @rdname BatchCore-methods
##' @keywords internal
setMethod(f = 'BatchCore',
          signature = c(p = 'big.matrix'),
          definition = function(p, idx, FUN, ..., n = 1) {

            ## register multiple core
            cl <- makeCluster(n)
            registerDoParallel(cl)

            ## describe files
            idxBigDesc <- describe(p)

            ## share package
            worker.init <- function (packages) {
              for (p in packages) {
                library(p, character.only=TRUE)
              }
              return(NULL)
            }
            clusterCall(cl, worker.init, 'bigmemory')

            ppiNum <- nrow(idx)

            ## get each value
            EachCall <- function(i,
                                 descFile,
                                 FUN,
                                 ...) {
              ## print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
              idxBig <- attach.big.matrix(descFile)

              f <- p[idxBig[i, 1], ]
              t <- p[idxBig[i, 2], ]
              eachVal <- FUN(eachArg = list(f = f, t = t, uniID = i), ...)

              gc()

              return(eachVal)
            }

            batchVec <- parSapply(cl,
                                  1:ppiNum,
                                  EachCall,
                                  idxBigDesc,
                                  FUN,
                                  ...)

            stopCluster(cl)

            return(batchVec)
          })
