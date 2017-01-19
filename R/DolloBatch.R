##' Dollo's parsimony distance of input phylogenetic profiles
##'
##' DolloDistBatch(): Dollo's parsimony distance in batch mode.
##'
##' @title Batch process of Dollo's parsimony distance.
##' @inheritParams SimDistBatch
##' @inheritParams DolloDist
##' @return A numeric vector.
##' @examples
##' require('ape')
##'
##' descPath <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfile")
##' tree <- read.nexus(descPath)
##' pathList <- nodepath(tree)
##'
##' pair1 <- system.file('extdata', 'bioinfoGenepair1.txt', package = "PhyloProfile")
##' oneP <- read.table(pair1, sep = '\t', row.names = 1)
##' sampleP <- matrix(sample(0:1, nrow(oneP) * 1000, replace = TRUE), ncol = 1000)
##' colnames(sampleP) <- paste0('gene', 1:1000)
##' rownames(sampleP) <- rownames(oneP)
##' ft <- matrix(paste0('gene', sample(1:1000, 10)), ncol = 2)
##' rownames(ft) <- paste0('p_', 1:nrow(ft))
##'
##' DolloDistBatch(ftMat = ft, profileMat = sampleP, edgeMat = tree$edge, tipPath = pathList, n = 2)
##' @seealso dollo
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @export
##' 
DolloDistBatch <- function(ftMat, profileMat, edgeMat, tipPath, n = 1) {

  ## register multiple core
  registerDoParallel(cores = n)

  ppiNames <- rownames(ftMat)
  ppiNum <- nrow(ftMat)

  geneNames <- colnames(profileMat)

  batchVec <- foreach(i = 1:ppiNum, .combine = c) %dopar% {
    print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
    genepair <- profileMat[, geneNames %in% ftMat[i, 1:2]]
    eachSD <- DolloDist(edgeMat, tipPath, genepair[, 1], genepair[, 2])

    return(eachSD)
  }

  ## stop multiple core
  stopImplicitCluster()

  return(batchVec)
}
