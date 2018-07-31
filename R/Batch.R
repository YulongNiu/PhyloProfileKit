##' @include AllClasses.R AllGenerics.R
NULL

##' Parallel framework for \code{PPIdx} and \code{PPTreeIdx}.
##'
##' \code{Batch()}: Parallel framework to analysis paired profiles in batch mode.
##'
##' @title Parallel framework
##' @inheritParams Batch
##' @return A \code{numeric vector}.
##' @examples
##' require('magrittr')
##' require('ape')
##'
##' tree <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfileKit") %>% read.nexus
##' ppPath <- system.file('extdata', 'bioinfoProfile.csv', package = "PhyloProfileKit")
##'
##' sceP <- ppPath %>% read.csv(row.names = 1) %>% as.matrix %>% PP
##' scePI <- PPIdx(sceP, 1:6, 1:6)
##' scePTI <- sceP %>% PPTree(tree) %>% PPTreeIdx(1:6, 1:6)
##'
##' ## Mutual information
##' Batch(scePI, method = 'SimCor', n = 2)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom RcppParallel setThreadOptions
##' @importFrom bigmemory is.big.matrix
##' @rdname Batch-methods
##' @exportMethod Batch
##'
setMethod(f = 'Batch',
          signature = c(x = 'PPIdx'),
          definition = function(x, method, ..., n) {

            ## check method
            ms <- c('SimCor', 'SimJaccard', 'SimMIBin', 'SimMIConti',
                   'DistHamming', 'DistManhattan', 'DistEuclidean',
                   'DistMinkowski', 'custom')
            midx <- pmatch(method, ms)

            if (is.na(midx)) {
              stop('Invalid similarity/distance method')
            } else {
              m <- ms[midx]
            }

            ## check arguments
            args <- list(...)
            if (method == 'custom') {
              funcPtr = arguments[["func"]]
              if (is.null(funcPtr)) {
                stop('Parameter "func" is missing.')
              } else {}
            } else {}

            ## set parallel threads
            setThreadOptions(numThreads = n)

            ## parallel idx
            p <- PPData(x)
            idx <- IdxData(x)

            if (is.big.matrix(idx)) {
              bv <- BatchBigmat(p, idx, list(method = m), args)
            } else {
              bv <- BatchMat(p, idx, list(method = m), args)
            }

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = m)

            return(bvRes)
          })


