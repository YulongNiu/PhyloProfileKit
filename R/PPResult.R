##' @include PP.R
NULL

##' @inheritParams show
##' @importFrom utils str
##' @importFrom methods show
##' @rdname show-methods
##' @exportMethod show
##
setMethod(f = 'show',
          signature = 'PPResult',
          definition = function(object){

            res <- object@.Data

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profiling results"\n')
            cat('class: ', class(object), '\n')
            cat('method: ', object@method, '\n')
            cat('#linkages: ', length(res), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(res)
            str(object@idx)
            str(object@pnames)
            })
