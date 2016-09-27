##' BayesTraits
##'
##' BayesTraitsMode(): Test the functinal linkages by BayesTraits using independent or dependent mode. 
##'
##' BayesTraits(): Test the a pair of functional linkage.
##'
##' BayesTraitsBatch(): Test a batch of functional linkages
##'
##' @title Use BayesTraits to test the functional linkages.
##' @param BayesTraitsPath Full path of the program BayesTraits.
##' @param treeFilePath Full path of tree file, and the the tree should be rooted and in nexus format.
##' @param binFilePath Full path of phylogenetic profiles, and only contains two genes. The speperate character by tab. If the paired genes are all "1" or all "0", then 0 will be returned.
##' @param uniID Unique ID for output file, and it is designed for parallel programming.
##' @param priorAll Prior distribution.
##' @param iterNum Iteration number only for 'MCMC' method, the default value is 1010000.
##' @param method Either 'MCMC' or 'ML' (Maximum Likelihood) method.
##' @param mode Either 'Indep' for independent and 'Dep' for dependent.
##' @return
##' 
##' BayesTraitsMode(): Harmonic mean.
##'
##' BayesTraits(): A list.
##'
##' BayesTraitsBatch(): A list.
##' 
##' @examples
##'
##' descPath <- system.file('extdata', package = "PhyloProTree")
##' 
##' tree <- file.path(descPath, 'bioinfoTree.nex')
##' pair00 <- file.path(descPath, 'bioinfoGenepair00.txt')
##' pair11 <- file.path(descPath, 'bioinfoGenepair11.txt')
##' pair1 <- file.path(descPath, 'bioinfoGenepair1.txt')
##'
##' BTtree <- file.path(descPath, 'Primates.trees')
##' BTpair <- file.path(descPath, 'Primates.txt')
##'
##' oneP <- read.table(pair1, sep = '\t', row.names = 1)
##' sampleP <- matrix(sample(0:1, nrow(oneP) * 1000, replace = TRUE), ncol = 1000)
##' colnames(sampleP) <- paste0('gene', 1:1000)
##' rownames(sampleP) <- rownames(oneP)
##' ft <- matrix(paste0('gene', sample(1:1000, 10)), ncol = 2)
##' rownames(ft) <- paste0('p_', 1:nrow(ft))
##' 
##' \dontrun{
##' ## replace "BTPath" with the full path of BayesTraits in your computer
##' ## all absence in the profile
##' BayesTraits(BTPath, tree, pair00)
##' ## all presence in the profile
##' BayesTraits(BTPath, tree, pair11)
##' ## a nomal profile
##' BayesTraits(BTPath, tree, pair1)
##'
##' ## examples from BayesTraits
##' BayesTraits(BTPath, BTtree, BTpair)
##' BayesTraits(BTPath, BTtree, BTpair, method = 'MCMC')
##'
##' ## batch
##' BayesTraitsBatch(ftMat = ft, profileMat = sampleP, n = 2,
##' BayesTraitsPath = BTPath, treeFilePath = BTtree)
##' }
##' 
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.evolution.reading.ac.uk/BayesTraits.html}
##' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/?term=17090580}
##' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/?term=16103904}
##' @importFrom utils read.table write.table
##' @rdname Bayes
##' @export
BayesTraitsMode <- function(BayesTraitsPath, treeFilePath, binFilePath, method = 'ML', priorAll = 'exp 10', iterNum = 1010000, mode = 'Indep', uniID = '1') {

  ## check if all present or all absent
  profile <- read.table(binFilePath, sep = '\t', row.names = 1)

  if (sum(profile) == 0 || sum(profile) == 2 * nrow(profile)) {
    return(0)
  } else {}

  ## creat order file
  if (mode == 'Indep') {
    modeType <- 2
  }
  else if (mode == 'Dep') {
    modeType <- 3
  }

  if (method == 'MCMC') {
    ## simulation type is 'MCMC'
    simuType <- '2'
    priorType <- paste('PriorAll', priorAll)
    iterNum <- paste('Iterations', iterNum)
    orderMat <- matrix(c(modeType, simuType, priorType, iterNum, 'Run'), ncol = 1)
  }
  else if (method == 'ML') {
    simuType <- '1'
    orderMat <- matrix(c(modeType, simuType, 'Run'), ncol = 1)
  }
  orderFileName <- paste('orderFile_', uniID, '.txt', sep = '')
  write.table(orderMat, row.names = FALSE, col.names = FALSE, quote = FALSE, file =orderFileName)

  ## run BayesTraits
  runFileName <- paste('run_', uniID, '.txt', sep = '')
  system(paste(BayesTraitsPath, treeFilePath, binFilePath, '<', orderFileName, '>', runFileName))

  ## get harmonic mean
  if (method == 'MCMC') {
    harMean <- system(paste("sed -e '$!{h;d;}' -e x ", runFileName, " | awk '{print $3}'", sep = ''), intern = TRUE)
  }
  else if (method == 'ML') {
    harMean <- system(paste("sed -e '$!{h;d;}' -e x ", runFileName," | awk '{print $2}'", sep = ''), intern = TRUE)
  }
  harMean <- as.numeric(harMean)

  ## rm file
  interFileName <- dir(pattern = paste(binFilePath, '.log*', sep = ''))
  moveFileName <- c(orderFileName, runFileName, interFileName)
  file.remove(moveFileName)
  
  return(harMean)
}



##' @param ... Inherit parameters from the "BayesTraitsMode()" in this package.
##' @importFrom stats pchisq
##' @rdname Bayes
##' @export
BayesTraits <- function(...){
  
  indepMean <- BayesTraitsMode(..., mode = 'Indep')
  depMean <- BayesTraitsMode(..., mode = 'Dep')

  ## chisqure test
  ## 2 * (Mean_dep - Mean_indep)
  logFactor <- 2 * (depMean - indepMean)
  
  ## not right
  ## pValue <- pchisq(logFactor, df = 4)
  pValue <- 1 - pchisq(logFactor, df = 4)

  bayesRes <- list(indepMean = indepMean,
                   depMean = depMean,
                   logFactor = logFactor,
                   p = pValue)
  
  return(bayesRes)
} 


##' @inheritParams SimDistBatch
##' @inheritParams BayesTraits
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom utils write.table
##' @rdname Bayes
##' @export
BayesTraitsBatch <- function(ftMat, profileMat, n = 1, ...) {

  ## register multiple core
  registerDoParallel(cores = n)

  ppiNames <- rownames(ftMat)
  ppiNum <- nrow(ftMat)

  geneNames <- colnames(profileMat)

  pathList <- foreach(i = 1:ppiNum) %dopar% {
    print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
    genepair <- profileMat[, geneNames %in% ftMat[i, 1:2]]

    ## write local file
    filename <- paste('genepair', i, '.txt', sep = '')
    write.table(genepair, col.names = FALSE, quote = FALSE, sep = '\t', file = filename)
    
    logLR <- BayesTraits(binFilePath = filename, uniID = i, ...)

    ## remove local file
    file.remove(filename)
    
    return(logLR)
  }

  names(pathList) <- ppiNames

  ## stop multiple core
  stopImplicitCluster()

  return(pathList)
}

## BayesTraits('/home/Yulong/Biotools/BayesTraitsV2/BayesTraitsV2',
##             '/home/Yulong/RESEARCH/NewRpkg/ppt/PhyloProTree/inst/extdata/bioinfoTree.nex',
##             '/home/Yulong/RESEARCH/NewRpkg/ppt/PhyloProTree/inst/extdata/bioinfoGenepair1.txt',
##             method = 'ML')

## oneP <- read.table('/home/Yulong/RESEARCH/NewRpkg/ppt/PhyloProTree/inst/extdata/bioinfoGenepair1.txt', sep = '\t', row.names = 1)
## sampleP <- matrix(sample(0:1, nrow(oneP) * 1000, replace = TRUE), ncol = 1000)
## colnames(sampleP) <- paste0('gene', 1:1000)
## rownames(sampleP) <- rownames(oneP)
## ft <- matrix(paste0('gene', sample(1:1000, 10)), ncol = 2)
## rownames(ft) <- paste0('p_', 1:nrow(ft))

## BayesTraitsBatch(ftMat = ft, profileMat = sampleP, n = 2, BayesTraitsPath = '/home/Yulong/Biotools/BayesTraitsV2/BayesTraitsV2', treeFilePath = '/home/Yulong/RESEARCH/NewRpkg/ppt/PhyloProTree/inst/extdata/bioinfoTree.nex')
