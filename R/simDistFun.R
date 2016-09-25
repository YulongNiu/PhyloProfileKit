##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param ftMat 
##' @param profileMat 
##' @param FUN 
##' @param n 
##' @return 
##' @examples 
##' @author Yulong Niu \email{niuylscu@@gmail.com}
SimDistBatch <- function(ftMat, profileMat, FUN, n = 1) {
  
  ## register multiple core
  registerDoParallel(cores = n)

  
  ## stop multiple core
  stopImplicitCluster()
}



##' Correlation or distance between a pair of phylogenetic profile
##'
##' SimCor(): Person's correlation coefficient.
##' SimJaccard(): Jaccard similarity
##' SimMI(): Mutual information
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
##' ## Mutual information
##' MIAB <- SimMI(ab)
##' 
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



##' @inheritParams SimCor
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname simdist
##' @export
##'
##' 
SimMI <- function(pairProfile) {
  
  combVec <- pairProfile[1, ] + 2 * pairProfile[2, ]
  
  N <- ncol(pairProfile)
  A <- sum(combVec == 3)
  B <- sum(combVec == 1)
  C <- sum(combVec == 2)
  D <- N - A - B - C

  eachMI <- function(p1, p2, p3, n) {
    eachI <- p1 * log(n * p1 / ((p1 + p2) * (p1 + p3))) / n
    return(eachI)
  }

  I <- eachMI(A, B, C, N) + eachMI(B, A, D, N) + eachMI(C, A, D, N) + eachMI(D, C, B, N)

  return(I)
}
