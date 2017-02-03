##' @include PPIdx.R
NULL

##' @inheritParams show
##' @importFrom utils str
##' @importFrom methods show
##' @importFrom ape Ntip
##' @rdname show-methods
##' @exportMethod show
##
setMethod(f = 'show',
          signature = 'PPTreeIdx',
          definition = function(object){

            p <- PPData(object)
            idx <- object@idx
            tree <- object@tree

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profile with linkage indices"\n')
            cat('class: ', class(object), '\n')
            if (isBinMat_internal(p)) {
              cat('profile: "binning"', '\n')
            } else {
              cat('profile: "continuous"', '\n')
            }
            cat('#species: ', ncol(p), '\n')
            cat('#proteins: ', nrow(p), '\n')
            cat('#linkages: ', nrow(idx), '\n')
            cat('#tips: ', Ntip(tree), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(object)

            })

##' The constructor the \code{PPTreeIdx} class
##'
##' Construct a \code{PP} object. The species (columns) are in the same order of phylogenetic tree tips, and for the species not in the tree are deleted.
##'
##' @title Constructor of \code{PPTreeIdx}
##' @param pidx A \code{PPIdx} object.
##' @param tree A \code{phylo} object.
##' @return A \code{PPTreeIdx} object.
##' @examples
##' require('magrittr')
##' require('ape')
##'
##' ppLink <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PP %>% PPIdx(1:5, 1:5)
##'
##' ppTree <- rtree(10, tip.label = paste0('spe', 10:1))
##' PPTreeIdx(ppLink, ppTree)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>% %<>% extract
##' @export
##'
PPTreeIdx <- function(pidx, tree) {

  p <- PPData(pidx)

  ## compare species with tree tips
  p %<>% colnames %>% match(tree$tip.label, .) %>% extract(p, , ., drop = FALSE)

  return(new('PPTreeIdx', p, idx = pidx@idx, tree = tree))
}
