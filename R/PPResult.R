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


##' Output of profiling results.
##'
##' Convert profiling results into a data frame with protein/gene names.
##'
##' @title Profiling results
##' @param x A \code{PPResult} object.
##' @param row.names 
##' @param optional
##' @param ...
##' @return data.frame of profiling results.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @seealso \code{\link[base]{as.data.frame}}
##' @exportMethod as.data.frame
##'
setMethod(f = 'as.data.frame',
          signature = c(x = 'PPResult'),
          definition = function(x, row.names = NULL, optional = FALSE, ...) {
            pn <- x@pnames
            idx <- x@idx
            res <- data.frame(pn[idx[, 1]],
                              pn[idx[, 2]],
                              x@.Data)
            rownames(res) <- rownames(idx)
            colnames(res) <- c(colnames(idx), x@method)

            return(res)
          })
