

##' Correlation or distance between a pair of phylogenetic profile
##'
##' SimCor(): Person's correlation coefficient.
##' SimJaccard(): Jaccard similarity
##'
##' 
##' @title Correlation and distance
##' @param pairProfile A paired phylogenetic profile. Names of rows are genes and names of columns are species
##' @return A numeric value.
##' @examples
##' ## alpha and beta subunits from the F-type ATP synthase.
##' data(fatp)
##' ab <- fatp$atpPhylo[c('ATP5A1', 'ATP5B'), ]
##'
##' ## Person's correlation coefficient
##' corAB <- SimCor(ab)
##' ## Jaccard similarity
##' jacAB <- SimJaccard(ab)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom stats cor
##' @rdname simdist
##' @export
##'
##' 
SimCor <- function(pairProfile) {
  return(cor(pairProfile[1, ], pairProfile[2, ]))
}



##' @inheritParams SimCor
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom vegan vegdist
##' @rdname simdist
##' @export
##'
##' 
SimJaccard <- function(pairProfile) {
  return(1 - vegdist(pairProfile, method = 'jaccard'))
}



SimMI <- function(pairProfile) {
  
}
