library('ape')
context('dollo')

InferGainNodesR <- function(gainList) {

  ## The latest common ancestor for all gain tips
  commonAncestorVec <- Reduce(intersect, gainList)
  commonAncestor <- commonAncestorVec[length(commonAncestorVec)]

  gainVec <- lapply(gainList, function(x) {
    eachGain <- x[which(x == commonAncestor) : length(x)]
    return(eachGain)
  })

  gainVec <- unique(unlist(gainVec))

  return(gainVec)
  
}

InferEdgeR <- function(tree, tipPath, pr) {

  ## initial edge matrix
  edges <- tree$edge
  edgeNum <- nrow(edges)

    ## no gains
  if (sum(pr) == 0) {
    edgeMat <- matrix(rep(0, 2 * edgeNum), ncol = 2)
    return(edgeMat)
  } else {}

  ## gain tips
  gainTips <- tipPath[which(pr == 1)]

  ## infer gain nodes and tips
  gains <- InferGainNodesR(gainTips)

  edgeMat <- cbind(as.numeric(edges[, 1] %in% gains),
                   as.numeric(edges[, 2] %in% gains))
  return(edgeMat)
}

DolloDistR <- function(tree, tipPath, pr1, pr2) {
  
  edgeMat <- cbind(InferEdgeR(tree, tipPath, pr1),
                   InferEdgeR(tree, tipPath, pr2))

  dist <- sum(abs((edgeMat[, 1] - edgeMat[, 2]) - (edgeMat[, 3] - edgeMat[, 4])))

  return(dist)
}


testTreeText <- '((((t1, t2),(t3, t4)),t5), (t6, (t7, t8)));'
testTree <- read.tree(text = testTreeText)
pathList <- nodepath(testTree)

#####################Two examples from the 17535793 paper###############
test_that('Dollo distance is 2.', {
  expect_equal(DolloDist(testTree$edge, pathList, c(0, 0, 0, 0, 1, 1, 1, 1), c(0, 0, 0, 0, 0, 1, 1, 1)), 2)
})

test_that('Dollo distance is 1.', {
  expect_equal(DolloDist(testTree$edge, pathList, c(1, 1, 1, 1, 1, 0, 0, 1), c(0, 0, 0, 0, 1, 0, 0, 1)), 1)
})
########################################################################

#####################test R version and Arma version####################
testNum <- 10000
gainMat1 <- matrix(sample(0:1, testNum * 8, replace = TRUE), ncol = testNum)
gainMat2 <- matrix(sample(0:1, testNum * 8, replace = TRUE), ncol = testNum)

dollo1 <- numeric(testNum)
dollo2 <- numeric(testNum)
for(i in 1:testNum){dollo1[i] <- DolloDistR(testTree, pathList, gainMat1[, i], gainMat2[, i])}
for(i in 1:testNum){dollo2[i] <- DolloDist(testTree$edge, pathList, gainMat1[, i], gainMat2[, i])}

test_that('Dollo distances are equal in two versions.', {
  expect_equal(sum(dollo1 == dollo2), testNum)
})
#######################################################################



