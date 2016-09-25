##' Similarity or distance of input phylogenetic profiles
##'
##' SimDistBatch(): Similarity and distance in batch mode.
##' 
##' @title Batch process of similarity and distance.
##' @param ftMat A two column matrix which should at least have rownames.
##' @param profileMat The phylogenetic profile data with 1 and 0 denoting the presence and absence of orthologous, respectively. It is a named numeric matrix, columns are species and rows are genes.
##' @param FUN Functions to calculate single similarity or distance.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @return A numeric vector
##' @examples
##' data(fatp)
##' f1 <- t(combn(rownames(fatp$atpPhylo)[1:6], 2))
##' SimDistBatch(f1, fatp$atpPhylo, SimCor, 2)
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

  geneNames <- rownames(profileMat)

  batchVec <- foreach(i = 1:ppiNum, .combine = c) %dopar% {
    print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
    genepair <- profileMat[geneNames %in% ftMat[i, 1:2], ]
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
  return(cor(pairProfile[1, ], pairProfile[2, ]))
}



##' @inheritParams SimCor
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname simdist
##' @export
##'
##' 
SimJaccard <- function(pairProfile) {
  
  f <- pairProfile[1, ]
  t <- pairProfile[2, ]

  A <- sum((f + 2*t) == 3)

  jac <- A / (sum(f) + sum(t) - A)
  
  return(jac)
}



##' @inheritParams SimCor
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname simdist
##' @export
##'
##' 
DistHamming <- function(pairProfile) {
  return(sum(pairProfile[1, ] != pairProfile[2, ]))
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

  NaN2Zero <- function(x) {
    if (is.na(x)) {
      x <- 0
    } else {}

    return(x)
  }

  I <- NaN2Zero(eachMI(A, B, C, N)) + NaN2Zero(eachMI(B, A, D, N)) + NaN2Zero(eachMI(C, A, D, N)) + NaN2Zero(eachMI(D, C, B, N))

  return(I)
}
