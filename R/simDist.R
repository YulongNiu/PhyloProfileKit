##' Similarity or distance of input phylogenetic profiles
##'
##' SimDistBatch(): Similarity and distance in batch mode.
##' 
##' @title Batch process of similarity and distance.
##' @param ftMat A two column matrix which should at least have rownames.
##' @param profileMat The phylogenetic profile data with 1 and 0 denoting the presence and absence of orthologous, respectively. It is a named numeric matrix, columns are genes and rows are species.
##' @param FUN Functions to calculate single similarity or distance.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @return A numeric vector
##' @examples
##' data(fatp)
##' f1 <- t(combn(rownames(fatp$atpPhylo)[1:6], 2))
##' SimDistBatch(f1, t(fatp$atpPhylo), SimCor, 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @seealso simdist
##' @export
##' 
SimDistBatch <- function(ftMat, profileMat, FUN, n = 1) {
  
  ## register multiple core
  registerDoParallel(cores = n)

  ppiNames <- rownames(ftMat)
  ppiNum <- nrow(ftMat)

  geneNames <- colnames(profileMat)

  batchVec <- foreach(i = 1:ppiNum, .combine = c) %dopar% {
    print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
    genepair <- profileMat[, geneNames %in% ftMat[i, 1:2]]
    eachSD <- FUN(genepair)

    return(eachSD)
  }
  
  ## stop multiple core
  stopImplicitCluster()

  return(batchVec)
}



##' Similarity or distance of paired phylogenetic profile
##'
##' SimCor(): Person's correlation coefficient.
##' SimJaccard(): Jaccard similarity.
##' SimMI(): Mutual information.
##' DistHamming(): Hamming distance.
##' 
##' @title similarity and distance
##' @param pairProfile A paired phylogenetic profile, columns are genes and rows are species.
##' @return A numeric value.
##' @examples
##' ## alpha and beta subunits from the F-type ATP synthase.
##' data(fatp)
##' ab <- t(fatp$atpPhylo[c('ATP5A1', 'ATP5B'), ])
##'
##' ## Person's correlation coefficient
##' corAB <- SimCor(ab)
##' ## Jaccard similarity
##' jacAB <- SimJaccard(ab)
##' ## Mutual information
##' MIAB <- SimMI(ab)
##' ## Hamming distance
##' hamAB <- DistHamming(ab)
##' 
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom stats cor
##' @rdname simdist
##' @seealso SimDistBatch
##' @export
##'
##' 
SimCor <- function(pairProfile) {
  return(cor(pairProfile[, 1], pairProfile[, 2]))
}

## library('Rcpp')
## library('RcppArmadillo')
## library('microbenchmark')
## library('bioDist')
## load('../data/fatp.RData')
## sourceCpp('../src/simDistCpp.cpp')

## microbenchmark(
##   'R' = for (i in 1:1000) {SimMIR(t(fatp$atpPhylo[sample(1:17, 2, replace = TRUE), ]))},
##   'arma' = for (i in 1:1000) {SimMI(t(fatp$atpPhylo[sample(1:17, 2, replace = TRUE), ]))},
##   'bioDist' = for (i in 1:1000) {mutualInfo(fatp$atpPhylo[sample(1:17, 2, replace = TRUE), ])}
## )
## for(i in 1:10) {
##   ab <- t(fatp$atpPhylo[sample(1:17, 2, replace = TRUE), ])
##   stopifnot(
##     all.equal(SimMIR(ab), as.numeric(mutualInfo(t(ab)))),
##     all.equal(SimMIR(ab), SimMI(ab))
##   )
## }
