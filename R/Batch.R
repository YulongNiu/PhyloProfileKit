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
##' ## Person's correlation coefficient
##' Batch(scePI, method = 'SimCor', n = 2)
##'
##' ## Jaccard similarity
##' Batch(scePI, method = 'SimJaccard', n = 2)
##'
##' ## Mutual information
##' Batch(scePI, method = 'SimMI', n = 2)
##'
##' ## Hamming distance
##' Batch(scePI, method = 'DistHamming', n = 2)
##'
##' ## Manhattan distance
##' Batch(scePI, method = 'DistManhattan', n = 2)
##'
##' ## Euclidean distance
##' Batch(scePI, method = 'DistEuclidean', n = 2)
##'
##' ## Minkowski distance
##' Batch(scePI, method = 'DistMinkowski', p = 4, n = 2)
##'
##' ## custom distance
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

            p <- PPData(x)
            idx <- IdxData(x)

            ## check method
            ms <- c('SimCor', 'SimJaccard', 'SimMI',
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
            }
            else if (method == 'SimMI') {
              m <- ifelse(isBinMat_(PP), 'SimMIBin', 'SimMIConti')
            } else {}

            ## set parallel threads
            setThreadOptions(numThreads = n)

            ## parallel idx
            if (is.big.matrix(idx)) {
              bv <- BatchBigmat(p, idx, list(method = m), args)
            } else {
              bv <- BatchMat(p, idx@address, list(method = m), args)
            }

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = m)

            return(bvRes)
          })


