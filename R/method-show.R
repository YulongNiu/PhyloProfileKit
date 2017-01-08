##' @include AllClasses.R
NULL


##' Show method for \code{PP} instances
##'
##' @title show methods
##' @param object A \code{PP} instance.
##' @return Show messages.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom utils str
##' @importFrom methods show
##' @rdname method-show
##' @exportMethod show
##'
setMethod(f = 'show',
          signature = 'PP',
          definition = function(object){
            d <- object@.Data

            ##~~~~~~~~~~~~~cat message~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profile"\n')
            cat('class: ', class(object), '\n')
            if (is.integer(object)) {
              cat('type: "binning"', '\n')
            } else {
              cat('type: "continuous"', '\n')
            }
            cat('#species: ', ncol(d), '\n')
            cat('#gene/protein: ', nrow(d), '\n')
            cat('---\n')
            str(d)
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          })

