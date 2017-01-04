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


##' @inheritParams SimCor
##' @rdname simdist
##' @importFrom bioDist mutualInfo
##' @export
SimMIConti <- function(pairProfile) {
  return(as.numeric(mutualInfo(t(pairProfile))))
}


## library('Rcpp')
## library('RcppArmadillo')
## library('microbenchmark')
## ## library('PhyloProfile')
## sourceCpp('../src/simDistCpp.cpp')

## microbenchmark(
##   'R' = for(i in 1:1000){SimCorR(matrix(rnorm(10000), ncol = 2))},
##   'arma2' = for(i in 1:1000){SimCor(matrix(rnorm(10000), ncol = 2))},
##   'arma' = for(i in 1:1000){SimCor2(matrix(rnorm(10000), ncol = 2))}
## )

