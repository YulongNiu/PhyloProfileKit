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
##' tree <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfile") %>% read.nexus
##' ppPath <- system.file('extdata', 'bioinfoProfile.csv', package = "PhyloProfile")
##'
##' sceP <- ppPath %>% read.csv(row.names = 1) %>% as.matrix %>% PP
##' sceT <- PPIdx(sceP, 1:6, 1:6) %>% PPTreeIdx(tree)
##' Dollo(sceT, n = 2)
##' 
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ape nodepath
##' @rdname Dollo-methods
##' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/?term=17535793}{Dollo's parsimony description}

##' @exportMethod Dollo
##'
setMethod(f = 'Dollo',
          signature = c(x = 'PPTreeIdx'),
          definition = function(x, ..., n = 1) {

            tree <- x@tree
            em <- tree$edge
            tp <- nodepath(tree)

            bv <- Batch(x = x,
                        FUN = DolloDist,
                        edgeMat = em,
                        tipPath = tp,
                        ...,
                        n = n)

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = 'Dollo')
          })
