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
            cat('description: "phylogenetic profile with tree and linkage indices"\n')
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
##' Construct a \code{PPTreeIdx} object. The species (columns) are in the same order of phylogenetic tree tips, and for the species not in the tree are deleted.
##'
##' @title Constructor of \code{PPTreeIdx}
##' @param pt A \code{PPTree} object.
##' @inheritParams PPIdx
##' @return A \code{PPTreeIdx} object.
##' @examples
##' require('magrittr')
##' require('ape')
##'
##' tree <- rtree(8, tip.label = paste0('spe', 8:1))
##' ppTree <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PPTree(tree)
##'
##' ## with self linkages
##' PPTreeIdx(ppTree, 1:3, 1:3, self = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>% %<>%
##' @seealso PPIdx
##' @export
##'
PPTreeIdx <- function(pt, x, ..., bigmat = FALSE) {

  pidx <- PPIdx(pt, x, ..., bigmat)

  p <- PPData(pt)
  tree <- pt@tree
  idx <- pidx@idx

  return(new('PPTreeIdx', p, tree = tree, idx = idx))
}


