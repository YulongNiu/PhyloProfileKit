##' @include utilities.R
NULL

##' This class represents the data structure of phylogenetic profile.
##'
##' @slot .Data An integer matrix or a numeric matrix, of which the rows are genes/proteins and columns are species. It validates the rownames and colnames of the matrix.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @exportClass PP
##' 
setClass(Class = 'PP',
         contains = 'matrix',
         validity = function(object) {
           d <- object@.Data
           valiMat_internal(d, 'PP')
         })

## ## validate nrow and ncol
## ## validate rownames
setClass(Class = 'Idx',
         contains = 'matrix',
         validity = function(object) {
           d <- object@.Data
           valiMat_internal(d, 'Idx')
         })


## Test <- function(x) {
##   if (!(is.numeric(x) &&
##         is.matrix(x))) {
##     return('haha')
##   } else {
##     TRUE
##   }
## }

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

