##' @include utilities.R
NULL

##' This class represents the data structure of phylogenetic profile.
##'
##' @slot .Data An integer matrix or a numeric matrix, of which the rows are genes/proteins and columns are species. It validates the rownames and colnames of the profile matrix.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @exportClass PP
##' 
setClass(Class = 'PP',
         contains = 'matrix',
         validity = function(object) {
           d <- object@.Data
           valiMat_internal(d, 'profile')
         })


setClass('big.matrix',
         slot = c(address = 'externalptr'))

setClassUnion(name = 'PPMat',
              member = c('matrix', 'big.matrix'))

##' This class represents the data structure of phylogenetic profile with linkage indices.
##'
##' @slot idx An integer matrix with two columns. It validates the rownames and colnames of the profile.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @exportClass PPIdx
##' 
setClass(Class = 'PPIdx',
         slots = c(idx = 'PPMat'),
         contains = 'PP',
         validity = function(object) {
           d <- object@idx
           valiMat_internal(d, 'indices')
         })


##' This class represents the data structure of phylogenetic profiling results.
##'
##' @slot .Data A numeric vector representing the profiling data.
##' @slot idx An integer matrix with two columns. It validates the rownames and colnames of the profile.
##' @slot pnames A character vector storing rownames (protein names) of the profile.
##' @slot method A character string indicating the profiling method.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @exportClass PPResult
##' 
setClass(Class = 'PPResult',
         slots = c(idx = 'PPMat', pnames = 'character', method = 'character'),
         contains = 'numeric')


structure(list(), class = 'phylo')

setOldClass('phylo')


##' This class represents the data structure of the phylogenetic profile with linkage indices and a phylogenetic tree.
##'
##' @slot tree A numeric vector representing the profiling data.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @seealso The raw definition of \code{\link[ape]{phylo}}.
##' @exportClass PPTreeIdx
##' 
setClass(Class = 'PPTreeIdx',
         slots = c(tree = 'phylo'),
         contains = 'PPIdx')
