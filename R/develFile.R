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
