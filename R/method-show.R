##' @include AllClasses.R
NULL


##' Show method for \code{PP} instances
##'
##' @title show methods
##' @param object A \code{PP} instance.
##' @return messages
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
            cat('#\n# phylogenetic profile\n#\n')
            cat('#...species number', '\t', ncol(d), '\n')
            cat('#...gene/protein number', '\t', nrow(d), '\n')
            if (is.integer(object)) {
              cat('#...type', '\t', 'binning', '\n')
            } else {
              cat('#...type', '\t', 'continuous', '\n')
            }
            str(d)
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          })

