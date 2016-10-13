## library('microbenchmark')
## library('igraph')
## library('ape')
## library('Rcpp')
## library('RcppArmadillo')
## library('ggtree')

## sourceCpp('../src/Dollo.cpp')
## source('PlotState.R')

## ##########################BioinfoTest##############################
## load('/home/Yulong/RESEARCH/PhyloProTreeTest/BioinfoData/yeastData.RData')

## ## tree
## bioinfoTree <- read.nexus('/home/Yulong/RESEARCH/PhyloProTreeTest/BioinfoData/tree_rooted.nex')
## pathList <- nodepath(bioinfoTree)

## ## rank profile
## profile <- yeastProfile[match(bioinfoTree$tip.label, rownames(yeastProfile)), ]

## fromID <- '6320111'
## toID <- '6321851'

## fromID <- '6323420'
## toID <- '6323886'

## fromID <- '37362655'
## toID <- '6319567'
## PlotState(bioinfoTree,
##           CombineNodeState(bioinfoTree,
##                            pathList,
##                            profile[, fromID],
##                            profile[, toID]))

## PlotState(bioinfoTree,
##           t(WeightEdge(bioinfoTree$edge,
##                        profile[, fromID],
##                        profile[, toID])))
## #################################################################

## ##############################review test##########################
## testTreeText <- '((((t1, t2),(t3, t4)),t5), (t6, (t7, t8)));'
## testTree <- read.tree(text = testTreeText)
## pathList <- nodepath(testTree)

## testNum <- 10000
## set.seed(123)
## gainMat1 <- matrix(sample(0:1, testNum * 8, replace = TRUE), ncol = testNum)
## set.seed(456)
## gainMat2 <- matrix(sample(0:1, testNum * 8, replace = TRUE), ncol = testNum)

## microbenchmark(
##   'dollo' = for(i in 1:testNum){DolloDist(testTree$edge, pathList, gainMat1[, i], gainMat2[, i])}
## )
## ######################################################################

## library('microbenchmark')
## library('ape')
## library('Rcpp')
## library('RcppArmadillo')

## sourceCpp('../src/Dollo.cpp')
## sourceCpp('../src/simDistCpp.cpp')

## set.seed(123456)
## testTree <- rtree(1000)
## pathList <- nodepath(testTree)
## testNum <- 1000

## set.seed(123123)
## gainMat1 <- matrix(sample(0:1, testNum * 1000, replace = TRUE), ncol = testNum)
## set.seed(456456)
## gainMat2 <- matrix(sample(0:1, testNum * 1000, replace = TRUE), ncol = testNum)

## microbenchmark(
##   'gain' = for(i in 1:testNum){DolloDist(testTree$edge, pathList, gainMat1[, i], gainMat2[, i])}
## )

