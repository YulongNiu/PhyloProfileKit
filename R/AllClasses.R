##' This class represents the data structure of phylogenetic profile.
##'
##' @slot .Data An integer matrix or a numeric matrix, of which the rows are genes/proteins and columns are species. It validates the rownames and colnames of the matrix.
##' @examples
##' ppBinning <- new('pp',
##'                  matrix(sample(0:1, 10 * 20, replace = TRUE),
##'                         ncol = 20,
##'                         dimnames = list(paste0('gene', 1:10),
##'                                         paste0('spe', 1:20))))
##' ppContinuous <- new('pp',
##'                     matrix(rnorm(10 * 20),
##'                            ncol = 20,
##'                            dimnames = list(paste0('gene', 1:10),
##'                                            paste0('spe', 1:20))))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom methods setClass
##' @exportClass pp
setClass(Class = 'pp',
         contains = 'matrix',
         validity = function(object) {
           d <- object@.Data
           if (sum(dim(d)) == 0) {
             return(TRUE)
           } else {
             if (is.null(rownames(d)) ||
                 is.null(rownames(d))) {
               warn <- 'The pp matrix (.Data slot) needs rownames or colnames.'
               return(warn)
             } else {return(TRUE)}
           }
         })

## ## validate nrow and ncol
## ## validate rownames
## setClass(Class = 'idx',
##          contains = 'matrix')

## tmp3 <- new('idx', matrix(sample(1:10, 3 * 2), ncol = 2))

## ## validate boundary
## ## validate integer
## setClass(Class = 'ppBin',
##          slots = c(ppData = 'pp', selectIdx = 'idx'))

## tmp4 <- new('ppBin', ppData = tmp1, selectIdx = tmp3)


## ## validate boundary
## ## validate numeric
## setClass(Class = 'ppCont',
##          slots = c(ppData = 'pp', selectIdx = 'idx'))
## tmp5 <- new('ppCont', ppData = tmp2, selectIdx = tmp3)

## setClassUnion(name = 'ppIdx',
##               members = c('ppBin', 'ppCont'))

