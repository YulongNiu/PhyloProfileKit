context('simdist')

library('stats')
library('bioDist')

SimCorR <- function(pairProfile) {
  return(cor(pairProfile[, 1], pairProfile[, 2]))
}

SimJaccardR <- function(pairProfile) {
  f <- pairProfile[, 1]
  t <- pairProfile[, 2]

  A <- sum((f + 2*t) == 3)

  jac <- A / (sum(f) + sum(t) - A)
  return(jac)
}


DistHammingR <- function(pairProfile) {
  return(sum(abs(pairProfile[, 1] - pairProfile[, 2])))
}


SimMIBinR <- function(pairProfile) {
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

DistEuclideanR <- function(pairProfile) {
  return(sqrt(sum((pairProfile[, 1] - pairProfile[, 2])^2)))
}

SimMIContiR <- function(pairProfile) {
  return(as.numeric(mutualInfo(t(pairProfile))))
}

#######################test R version and Arma version############
testNum <- 10000

sd1 <- numeric(testNum)
sd2 <- numeric(testNum)
sd3 <- numeric(testNum)

##~~~~~~~~~~~~~~~~~~~~~~~~~cor~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(rnorm(2 * speNum), ncol = 2)

  sd1[i] <- SimCorR(pptmp)
  sd2[i] <- SimCor(pptmp[, 1], pptmp[, 2])
}

test_that('Pearson correlation coefficient are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~continuous MI~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(rnorm(2 * speNum), ncol = 2)

  sd1[i] <- SimMIContiR(pptmp)
  sd2[i] <- SimMIConti(pptmp[, 1], pptmp[, 2], bin = 10)

}

test_that('Continuous MI similarity are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~binning MI~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(sample(0:1, 2 * speNum, replace = TRUE), ncol = 2)

  sd1[i] <- SimMIBinR(pptmp)
  sd2[i] <- SimMIBin(pptmp[, 1], pptmp[, 2])
}

test_that('Binning MI similarity are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Jaccard~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(sample(0:1, 2 * speNum, replace = TRUE), ncol = 2)

  sd1[i] <- SimJaccardR(pptmp)
  sd2[i] <- SimJaccard(pptmp[, 1], pptmp[, 2])
}

test_that('Jaccard similarity are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Hamming~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(rnorm(2 * speNum), ncol = 2)

  sd1[i] <- DistHammingR(pptmp)
  sd2[i] <- DistHamming(pptmp[, 1], pptmp[, 2])
}

test_that('Hamming distances are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Euclidean MI~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {
  speNum <- 20
  pptmp <- matrix(rnorm(2 * speNum), ncol = 2)

  sd1[i] <- DistEuclideanR(pptmp)
  sd2[i] <- DistEuclidean(pptmp[, 1], pptmp[, 2])
  sd3[i] <- as.numeric(dist(t(pptmp), method = 'euclidean'))
}

test_that('Euclidean distances are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
  expect_equal(all.equal(sd1, sd3), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################


