##' @include PP.R
NULL

##' @inheritParams show
##' @importFrom utils str
##' @importFrom methods show
##' @importFrom ape Ntip
##' @rdname show-methods
##' @exportMethod show
##
setMethod(f = 'show',
          signature = 'PPTree',
          definition = function(object){

            p <- PPData(object)
            tree <- object@tree

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profile with tree"\n')
            cat('class: ', class(object), '\n')
            if (isBinMat_internal(p)) {
              cat('profile: "binning"', '\n')
            } else {
              cat('profile: "continuous"', '\n')
            }
            cat('#species: ', ncol(p), '\n')
            cat('#proteins: ', nrow(p), '\n')
            cat('#tips: ', Ntip(tree), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(object)

            })

##' The constructor the \code{PPTree} class
##'
##' Construct a \code{PPTree} object from a numeric matrix.
##'
##' @title Constructor of \code{PPTree}
##' @param tree A \code{phylo} object.
##' @inheritParams PP
##' @return A \code{PPTree} object.
##' @examples
##' require('magrittr')
##' require('ape')
##'
##' ppMat <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20)
##' ppTree <- rtree(8, tip.label = paste0('spe', 8:1))
##' PPTree(pp, ppTree)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>% %<>%
##' @export
##' 
PPTree <- function(x, tree) {

  ## check pp and names
  p <- x %>% PP %>% PPData

  ## remove species not in the tree
  p %<>% colnames %>% match(tree$tip.label, .) %>% extract(p, , ., drop = FALSE)

  return(new('PPTree', p, tree = tree))
}
