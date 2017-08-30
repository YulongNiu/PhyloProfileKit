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

            p <- object@.Data
            idx <- object@idx
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
