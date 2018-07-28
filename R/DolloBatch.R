##' @include AllClasses.R AllGenerics.R Batch.R
NULL

##' Dollo's parsimony distance
##'
##' Dollo's parsimony distance of paired profiles.
##'
##' @inheritParams Dollo
##' @title Batch process of Dollo's parsimony distance
##' @return A \code{PPResult} object.
##' @examples
##' require('magrittr')
##' require('ape')
##'
##' tree <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfileKit") %>% read.nexus
##' ppPath <- system.file('extdata', 'bioinfoProfile.csv', package = "PhyloProfileKit")
##'
##' sceP <- ppPath %>% read.csv(row.names = 1) %>% as.matrix %>% PP
##' scePTI <- sceP %>% PPTree(tree) %>% PPTreeIdx(1:6, 1:6)
##' Dollo(scePTI, n = 2)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom ape nodepath
##' @rdname Dollo-methods
##' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/?term=17535793}{Dollo's parsimony description}
##' @exportMethod Dollo
##'
setMethod(f = 'Dollo',
          signature = c(x = 'PPTreeIdx'),
          definition = function(x, ..., n) {

            dollo_internal <- function(eachArg, edgeMat, tipPath, ...) {
              return(DolloDist(edgeMat = edgeMat,
                               tipPath = tipPath,
                               f = eachArg$f,
                               t = eachArg$t))
            }

            tree <- x@tree
            em <- tree$edge
            tp <- nodepath(tree)

            bv <- Batch(x = x,
                        FUN = dollo_internal,
                        edgeMat = em,
                        tipPath = tp,
                        ...,
                        n = n)

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = 'Dollo')

            return(bvRes)
          })
