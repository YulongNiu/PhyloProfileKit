##' @include AllClasses.R AllGenerics.R Batch.R
NULL


##' BayesTraits
##'
##' A wrapper of the BayesTraits test of paired profiles. Please download and install \href{http://www.evolution.reading.ac.uk/BayesTraits.html}{BayesTraits}.
##'
##' @inheritParams BayesTraits
##' @title Batch process of BayesTraits
##' @return A \code{PPResult} object.
##' @examples
##' require('magrittr')
##' require('ape')
##'
##' tree <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfile") %>% read.nexus
##' ppPath <- system.file('extdata', 'bioinfoProfile.csv', package = "PhyloProfile")
##'
##' sceP <- ppPath %>% read.csv(row.names = 1) %>% as.matrix %>% PP
##' sceT <- PPIdx(sceP, 1:6, 1:6) %>% PPTreeIdx(tree)
##'
##' \dontrun{
##' ## replace "BTPath" with the full path of BayesTraits in your system, for example in Linux/OS
##' ## BTPath <- '/path/to/BayesTraitsV2/BayesTraitsV2'
##' ## or in Windows
##' ## BTPath <- '/path/to/BayesTraitsV2/BayesTraitsV2.exe'
##' BayesTraits(sceT, BayesTraitsPath = BTPath, n = 2)
##' }
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ape write.nexus
##' @rdname BayesTraits-methods
##' @references \href{http://www.evolution.reading.ac.uk/BayesTraits.html}{BayesTraits program}
##' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/?term=17090580}
##' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/?term=16103904}
##' @exportMethod BayesTraits
##'
setMethod(f = 'BayesTraits',
          signature = c(x = 'PPTreeIdx'),
          definition = function(x, ..., n = 1) {

            bt_internal <- function(eachArg, treeFilePath, ...) {
              genepair <- cbind(eachArg$f, eachArg$t)

              ## write local file
              filename <- paste('genepair', eachArg$uniID, '.txt', sep = '')
              write.table(genepair, col.names = FALSE, quote = FALSE, sep = '\t', file = filename)

              logLR <- BayesTraitsTest(treeFilePath = treeFilePath,
                                       binFilePath = filename,
                                       uniID = eachArg$uniID,
                                       ...)$logFactor

              ## remove local file
              file.remove(filename)

              return(logLR)
            }

            ## write tree
            treeFile <- 'tree.nex'
            write.nexus(x@tree, file = treeFile)

            bv <- Batch(x = x,
                        FUN = bt_internal,
                        treeFilePath = treeFile,
                        ...,
                        n = n)

            ## remove tree file
            file.remove(treeFile)

            bvRes <- new('PPResult',
                         bv,
                         idx = x@idx,
                         pnames = rownames(x@.Data),
                         method = 'BayesTraits')

            return(bvRes)
          })


##' BayesTraits Mode
##'
##' \code{BayesTraitsMode()}: Test the functinal linkages by BayesTraits using independent or dependent mode.
##'
##' \code{BayesTraitsTest()}: Test the a pair of functional linkage.
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
##' \code{BayesTraitsMode()}: Harmonic mean.
##'
##' \code{BayesTraitsTest()}: A list.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom utils read.table write.table
##' @rdname Bayes
##' @keywords internal
##'
BayesTraitsMode <- function(BayesTraitsPath,
                            treeFilePath,
                            binFilePath,
                            method = 'ML',
                            priorAll = 'exp 10',
                            iterNum = 1010000,
                            mode = 'Indep',
                            uniID = '1') {

  ## check if all present or all absent
  profile <- read.table(binFilePath, sep = '\t', row.names = 1)

  if (all(profile == 1) ||
      all(profile == 0)) {
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
  write.table(orderMat, row.names = FALSE, col.names = FALSE, quote = FALSE, file = orderFileName)

  ## run BayesTraits
  runFileName <- paste('run_', uniID, '.txt', sep = '')
  system(paste(BayesTraitsPath, treeFilePath, binFilePath, '<', orderFileName, '>', runFileName))

  ## get harmonic mean
  resStr <- readLines(runFileName)
  dataVec <- unlist(strsplit(resStr[length(resStr) - 1], split = '\t', fixed = TRUE))
  if (method == 'MCMC') {
    harMean <- as.numeric(dataVec[3])
  }
  else if (method == 'ML') {
    harMean <- as.numeric(dataVec[2])
  }

  ## rm file
  interFileName <- dir(pattern = paste(binFilePath, '.log*', sep = ''))
  moveFileName <- c(orderFileName, runFileName, interFileName)
  file.remove(moveFileName)

  return(harMean)
}



##' @param ... Inherit parameters from the "BayesTraitsMode()" in this package.
##' @importFrom stats pchisq
##' @rdname Bayes
##' @keywords internal
##'
BayesTraitsTest <- function(...){

  indepMean <- BayesTraitsMode(..., mode = 'Indep')
  depMean <- BayesTraitsMode(..., mode = 'Dep')

  ## chisqure test
  ## 2 * (Mean_dep - Mean_indep)
  logFactor <- 2 * (depMean - indepMean)
  pValue <- 1 - pchisq(logFactor, df = 4)

  bayesRes <- list(indepMean = indepMean,
                   depMean = depMean,
                   logFactor = logFactor,
                   p = pValue)

  return(bayesRes)
}
