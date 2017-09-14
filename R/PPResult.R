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

            ##~~~~~~~~~~~~~head~~~~~~~~~~~~
            cat('---\n')
            cat('description: "phylogenetic profiling results"\n')
            cat('class: ', class(object), '\n')
            cat('method: ', object@method, '\n')
            cat('#linkages: ', length(object@.Data), '\n')
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(object)
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
##' @importFrom bigmemory head
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
head.PPResult <- function(x, n = 6L, ...) {

  headx <- new('PPResult',
               head(x@.Data, n, ...),
               idx = head(x@idx, n, ...),
               pnames = x@pnames,
               method = x@method)
  headx %>% as.data.frame %>% return
}


##' @method tail PPResult
##' @importFrom magrittr %>%
##' @importFrom bigmemory tail
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
tail.PPResult <- function(x, n = 6L, ...) {

  tailx <- new('PPResult',
               tail(x@.Data, n, ...),
               idx = tail(x@idx, n, ...),
               pnames = x@pnames,
               method = x@method)
  tailx %>% as.data.frame %>% return
}


##' Compare \code{PPResult} objects with a numeric value.
##'
##' @title Compare methods
##' @param e1 A \code{PPResult} object.
##' @param e2 A numeric value.
##' @return A \code{logic vector}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>%
##' @importFrom methods Compare callGeneric
##' @rdname Compare-methods
##' @exportMethod Compare
##
setMethod('Compare',
          signature = c(e1 = 'PPResult', e2 = 'numeric'),
          function(e1, e2) {
            callGeneric(e1@.Data, e2) %>% return
          })



##' @importFrom methods validObject
##' @importFrom magrittr %>% %T>%
##' @rdname select-methods
##' @exportMethod [
##' 
setMethod(f = '[',
          signature = c(x = 'PPResult', j = 'missing'),
          definition = function(x, i, j, ..., drop = FALSE) {
            x@.Data <- x@.Data[i]
            x@idx <- x@idx[i, , drop = FALSE]
            x %T>% validObject %>% return
          })

