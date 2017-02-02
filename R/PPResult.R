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


##' @method as.data.frame PPResult
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
as.data.frame.PPResult <- function(x, ...) {
  pn <- x@pnames
  idx <- x@idx
  res <- data.frame(pn[idx[, 1]],
                    pn[idx[, 2]],
                    x@.Data)
  rownames(res) <- rownames(idx)
  colnames(res) <- c(colnames(idx), x@method)

  return(res)
}

##' @method head PPResult
##' @importFrom magrittr %>%
##' @importFrom utils head
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
head.PPResult <- function(x, n = 6L, ...) {

  headx <- new('PPResult',
               head(x@.Data, n, ...),
               idx = head(x@idx, n, ...),
               pnames = head(x@pnames, n, ...),
               method = x@method)
  headx %>% as.data.frame %>% return
}


##' @method tail PPResult
##' @importFrom magrittr %>%
##' @importFrom utils tail
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
tail.PPResult <- function(x, n = 6L, ...) {

  tailx <- new('PPResult',
               tail(x@.Data, n, ...),
               idx = tail(x@idx, n, ...),
               pnames = tail(x@pnames, n, ...),
               method = x@method)
  tailx %>% as.data.frame %>% return
}
