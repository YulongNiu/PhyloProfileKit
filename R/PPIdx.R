##' @include PP.R
NULL

##' @inheritParams show
##' @importFrom utils str
##' @importFrom methods show
##' @rdname show-methods
##' @exportMethod show
##
setMethod(f = 'show',
          signature = 'PPIdx',
          definition = function(object){

            p <- object@.Data
            idx <- object@idx

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
            cat('---\n')
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            str(p)
            str(idx)

            })


##' The constructor the \code{PPIdx} class
##'
##' Construct a \code{PP} object from a character matrix.
##'
##' @title Constructor of \code{PPIdx}
##' @param value A character matrix with two columns. The linkages containing proteins not in \code{p} (profile) are removed.
##' @param p A \code{PP} object
##' @return A \code{PPIdx} object
##' @examples
##' require('magrittr')
##'
##' ## construct a PP object without dimnames
##' ppBinning <- sample(0:1, 10 * 20, replace = TRUE) %>% matrix(ncol = 20) %>% PP
##' linkM <- sample(1:30, 20 * 2, replace = TRUE) %>% paste0('protein', .) %>% matrix(ncol = 2)
##' ppBinIdx <- PPIdx(linkM, ppBinning)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom magrittr %>%
##' @export
##' 
PPIdx <- function(value, p) {

  ## whole proteins
  wp <- rownames(p@.Data)

  if (is.matrix(value) &&
      is.character(value)) {

    colSize <- ncol(value)
    rowSize <- nrow(value)

    ## check 0 row or 0 columns
    if (colSize == 0 ||
        rowSize == 0) {
      return(new('PPIdx', p, idx = value))
    } else {}

    ## check colnames
    if (is.null(colnames(value))) {
      colnames(value) <- c('From', 'To')
    } else {}

    ## check rownames
    if (is.null(rownames(value))) {
      rownames(value) <- paste0('link', seq_len(rowSize))
    } else{}

    ## check linkage proteins in wp, NA if not in wp
    vcheck <- c(value) %>% match(wp) %>% matrix(ncol = 2, dimnames = dimnames(value))
    hasLogic <- !(is.na(vcheck[, 1]) | is.na(vcheck[, 2]))
    vcheck <- vcheck[hasLogic, , drop = FALSE]

    return(new('PPIdx', p, idx = vcheck))
  } else {
    return(value)
  }
}

