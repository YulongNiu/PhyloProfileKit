##' @include AllClasses.R AllGenerics.R utilities.R Batch.R
NULL

##' Similarity or distance
##'
##' Similarity and distance of paired profiles. If the input is a \code{PPTreeIdx} object, the paired profile is collapsed according to the phylogenetic tree.
##'
##' @inheritParams SimDist
##' @title Batch process of similarity and distance
##' @return A \code{PPResult} object.
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
##' SimDist(scePI, 'SimMI', n = 2)
##' SimDist(scePTI, 'SimMI', n = 2)
##' 
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname SimDist-methods
##' @exportMethod SimDist
##'
setMethod(f = 'SimDist',
          signature = c(x = 'PPIdx'),
          definition = function(x, method, ..., n) {

            sd_internal <- function(eachArg, M, ...) {
              return(M(eachArg$f, eachArg$t, ...))
            }

            M <- ChooseSimDistFun(x, method)

            bv <- Batch(x = x,
                        FUN = sd_internal,
                        M = M,
                        ...,
                        n = n)

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = method)

            return(bvRes)
          })



##' @inheritParams SimDist
##' @rdname SimDist-methods
##' @exportMethod SimDist
##' 
setMethod(f = 'SimDist',
          signature = c(x = 'PPTreeIdx'),
          definition = function(x, method, ..., n) {

            ct_internal <- function(eachArg, edgeMat, tipNum, M, ...) {

              ftMat <- CollapseTree(edgeMat = edgeMat,
                                    tipNum = tipNum,
                                    f = eachArg$f,
                                    t = eachArg$t)
              fnew <- ftMat[, 3]
              tnew <- ftMat[, 4]

              return(M(fnew, tnew, ...))
            }

            tree <- x@tree
            em <- tree$edge
            tn <- Ntip(tree)
            M <- ChooseSimDistFun(x, method)

            bv <- Batch(x = x,
                        FUN = ct_internal,
                        edgeMat = em,
                        tipNum = tn,
                        M = M,
                        ...,
                        n = n)

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = method)

            return(bvRes)
          })



##' Similarity or distance function.
##'
##' @inheritParams ChooseSimDistFun
##' @title Choose SimDist function
##' @return A function.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname ChooseSimDistFun-methods
##' @keywords internal
##'
setMethod(f = 'ChooseSimDistFun',
          signature = c(x = 'PP'),
          definition = function(x, method, ...) {

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

            return(m)
          })

