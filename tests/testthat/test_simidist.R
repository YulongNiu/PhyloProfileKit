context('simdist')

library('stats')
library('bioDist')
library('magrittr')

BatchMatR <- function(p, pidx, FUN) {

  ppiNum <- nrow(pidx)
  res <- double(ppiNum)

  for (i in seq_len(ppiNum)) {
    ft <- p[c(pidx[i, 1], pidx[i, 2]), ]
    res[i] <- FUN(t(ft))
  }

  return(res)
}

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

DistManhattanR <- function(pairProfile) {
  return(sum(abs(pairProfile[, 1] - pairProfile[, 2])))
}

DistHammingR <- function(pairProfile) {
  return(sum(pairProfile[, 1] != pairProfile[, 2]))
}

DistEuclideanR <- function(pairProfile) {
    return(sqrt(sum((pairProfile[, 1] - pairProfile[, 2])^2)))
}

SimMIContiR <- function(pairProfile) {
  return(as.numeric(mutualInfo(t(pairProfile))))
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


#######################test R version and Arma version############
testNum <- 2
prow  <- 1e5
pcol  <- 1e3
idxnum <- 1e3

sd1 <- vector('list', testNum)
sd2 <- vector('list', testNum)

##~~~~~~~~~~~~~~~~~~~~~~~~~cor~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, SimCorR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'SimCor'), list())
}

test_that('Pearson correlation coefficient are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Jaccard~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, SimJaccardR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'SimJaccard'), list())
}


test_that('Jaccard similarity are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Hamming~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, DistHammingR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'DistHamming'), list())
}

test_that('Hamming distances are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Manhattan~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, DistManhattanR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'DistManhattan'), list())
}

test_that('Manhattan distances are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Euclidean~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, DistEuclideanR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'DistEuclidean'), list())
}

test_that('Euclidean distances are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~binning MI~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, SimMIBinR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'SimMIBin'), list())
}


test_that('Binning MI similarity are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~continuous MI~~~~~~~~~~~~~~~~~~~
for(i in 1:testNum) {

  p <- sample(0:1, prow * pcol, replace = TRUE) %>%
    matrix(nrow = prow)
  pidx <- sample(1:prow, idxnum * 2, replace = TRUE) %>%
    matrix(nrow = idxnum)

  sd1[[i]] <- BatchMatR(p, pidx, SimMIContiR)
  sd2[[i]] <- BatchMat(p, pidx, list(method = 'SimMIConti'), list(bin = 10))
}

test_that('Continuous MI similarity are equal in two versions.', {
  expect_equal(all.equal(sd1, sd2), TRUE)
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################


