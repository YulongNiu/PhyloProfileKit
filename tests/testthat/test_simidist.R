context('simdist')

## SimJaccardR <- function(pairProfile) {
##   return(1 - vegdist(pairProfile, method = 'jaccard'))
## }

SimMIR <- function(pairProfile) {
  
  combVec <- pairProfile[, 1] + 2 * pairProfile[, 2]
  
  N <- nrow(pairProfile)
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

#######################test R version and Arma version of MI############
testNum <- 10000

MI1 <- numeric(testNum)
MI2 <- numeric(testNum)

for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(sample(0:1, 2 * speNum, replace = TRUE), ncol = 2, nrow = speNum)

  MI1[i] <- SimMIR(pptmp)
  MI2[i] <- SimMI(pptmp)
}

test_that('MI distances are equal in two versions.', {
  expect_equal(sum(MI1 == MI2), testNum)
})
########################################################################


