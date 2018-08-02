##' ######################test plot#######################
##'   require('grid')
##'   require('ggplot2')
##'   require('reshape2')
##'   require('ggdendro')
##'   require('gridExtra')
##'   require('PhyloProfile')
##' 
##'   data(fatp)
##' 
##'   phyloData <- fatp$atpPhylo
##'   speCol = fatp$specol
##'   geneCol = fatp$geneco
##' 
##'   geneNameSize = 3
##'   geneNameCol = 'grey55'
##'   geneBetweenBlockCol = NA
##'   presentCol = 'steelblue'
##'   absentCol = 'grey91'
##' 
##'   ## cluster genes and species
##'   hcGene <- hclust(dist(phyloData), method = 'average')
##'   rowInd <- hcGene$order
##'   hcSpe <- hclust(dist(t(phyloData)), method = 'average')
##'   colInd <- hcSpe$order
##' 
##'   ## order 'phyloData'
##'   orderedPhyloData <- phyloData[rowInd,colInd]
##'   orderedColNames <- colnames(orderedPhyloData)
##'   orderedRowNames <- rownames(orderedPhyloData)
##'   breaksRow <- 1:length(orderedRowNames)
##' 
##'   ## order 'geneCol'
##'   orderedGeneCol <- geneCol[match(rownames(orderedPhyloData), names(geneCol))]
##' 
##'   ## order 'speCol'
##'   orderedSpeCol <- speCol[match(colnames(orderedPhyloData), names(speCol))]
##' 
##'   ## melt data for ggplot2
##'   colnames(orderedPhyloData) <- 1:ncol(orderedPhyloData)
##'   rownames(orderedPhyloData) <- 1:nrow(orderedPhyloData)
##'   orderedPhyloData <- melt(orderedPhyloData)
##'   orderedPhyloData <- data.frame(geneNames = orderedPhyloData[, 1], speNames = orderedPhyloData[, 2], apData = factor(orderedPhyloData[, 3]))


## library('Rcpp')
## library('magrittr')
## library('RcppParallel')
## library('bigmemory')
## library('ape')
## library('RcppXPtrUtils')
## sourceCpp('src/batch.cpp')
## source('R/AllClasses.R')
## source('R/AllGenerics.R')
## source('R/Batch.R')
## source('R/PP.R')
## source('R/PPIdx.R')
## source('R/utilities.R')
## source('R/PPResult.R')

## tree <- read.nexus('inst/extdata/bioinfoTree.nex')
## sceP <- read.csv('inst/extdata/bioinfoProfile.csv', row.names = 1) %>% as.matrix %>% PP
## scePI <- PPIdx(sceP, 1:6, 1:6)

## BatchMat(PPData(scePI), IdxData(scePI), list(method = 'SimCor'), list())
## BatchMat(PPData(scePI), IdxData(scePI), list(method = 'SimCorCollapse'), list(edgeMat=tree$edge, tipNum = Ntip(tree)))


## Batch(scePI, method = 'SimCor', n = 2)
## euclideanFuncPtr <- cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) { return sqrt(arma::accu(arma::square(A - B))); }", depends = c("RcppArmadillo"))
## all.equal(Batch(scePI, method = 'SDCustom', func = euclideanFuncPtr, n = 2)@.Data, Batch(scePI, method = 'DistEuclidean', n = 2)@.Data)

## all.equal(Batch(scePI, method = 'SDCustom', func = euclideanFuncPtr, edgeMat=tree$edge, tipNum = Ntip(tree), n = 2)@.Data, Batch(scePI, method = 'DistEuclidean', edgeMat=tree$edge, tipNum = Ntip(tree), n = 2)@.Data)


