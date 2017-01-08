##' This class represents the data structure of phylogenetic profile.
##'
##' @slot .Data An integer matrix or a numeric matrix, of which the rows are genes/proteins and columns are species. It validates the rownames and colnames of the matrix.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom methods setClass
##' @exportClass PP
setClass(Class = 'PP',
         contains = 'matrix',
         validity = function(object) {
           d <- object@.Data
           if (ncol(d) == 0 ||
               nrow(d) == 0){
             return(TRUE)
           } else {
             if (is.null(rownames(d)) ||
                 is.null(rownames(d))) {
               warn <- 'The PP matrix (.Data slot) needs rownames or colnames.'
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

