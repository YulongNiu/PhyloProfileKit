##' Plot phylogenetic profiles with gene and species clusters.
##'
##' A combination plot of phylogenetic profiles with genes and species clusters.
##' 
##' @title Plot phylogenetic profiles
##' @param phyloData The phylogenetic profile data with 1 and 0 denoting the presence and absence of orthologous, respectively. The "phyloData" should be a numeric matrix, of which the row is gene and column is species.The "phyloData"has row names and column names which will be used for the dendrogram of row and column.
##' @param geneNameSize The size of gene names label, and the default value is 3.
##' @param geneNameCol The colour of gene names, and the default value is "grey55".
##' @param geneBetweenBlockCol The space color between gene blocks, and the default value is "NA" meaning no space color. If the number of genes is samll, for example less than 20, setting it as 'white' is fine. 
##' @param presentCol The color of present 1, the default value is "steelblue".
##' @param absentCol The color of present 0, the default value is "grey91".
##' @param speCol A vector of colors with names of species, which are the same as colnames of "phyloData" (may not in the same order). 
##' @param geneCol A vector of colors with names of species, which are the same as rownames of "phyloData" (may not in the same order).
##' @param widthsShinkage The shinkage width vector.
##' @param heightsShinkage The shinkage width vector. 
##' @return A plot object.
##' @examples
##' data(fatp)
##' ATPphyloPlot <- PlotPhyloProfile(fatp$atpPhylo, speCol = fatp$specol, geneCol = fatp$genecol)
##' \dontrun{
##' pdf('FATPprofilePlot.pdf')
##' ATPphyloPlot <- PlotPhyloProfile(fatp$atpPhylo, speCol = fatp$specol, geneCol = fatp$genecol)
##' dev.off()
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_text geom_tile geom_segment geom_point scale_fill_manual labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string coord_flip
##' @importFrom ggdendro dendro_data.hclust segment
##' @importFrom gridExtra grid.arrange
##' @importFrom reshape2 melt
##' @export
##' 
PlotPhyloProfile <- function(phyloData,
                            geneNameSize = 3,
                            geneNameCol = 'grey55',
                            geneBetweenBlockCol = NA,
                            presentCol = 'steelblue',
                            absentCol = 'grey91',
                            speCol,
                            geneCol,
                            widthsShinkage = c(0.9, 0.7, 0.3, 7),
                            heightsShinkage = c(0.3, 7)) {

  ## require('grid')
  ## require('ggplot2')
  ## require('reshape2')
  ## require('ggdendro')
  ## require('gridExtra')

  ## cluster genes and species
  hcGene <- hclust(dist(phyloData), method = 'average')
  rowInd <- hcGene$order
  hcSpe <- hclust(dist(t(phyloData)), method = 'average')
  colInd <- hcSpe$order
  
  ## order 'phyloData'
  orderedPhyloData <- phyloData[rowInd,colInd]
  orderedColNames <- colnames(orderedPhyloData)
  orderedRowNames <- rownames(orderedPhyloData)
  breaksRow <- 1:length(orderedRowNames)

  ## order 'geneCol'
  orderedGeneCol <- geneCol[match(rownames(orderedPhyloData), names(geneCol))]
  
  ## order 'speCol'
  orderedSpeCol <- speCol[match(colnames(orderedPhyloData), names(speCol))]
  
  ## melt data for ggplot2
  colnames(orderedPhyloData) <- 1:ncol(orderedPhyloData)
  rownames(orderedPhyloData) <- 1:nrow(orderedPhyloData)
  orderedPhyloData <- melt(orderedPhyloData)
  orderedPhyloData <- data.frame(geneNames = orderedPhyloData[, 1], speNames = orderedPhyloData[, 2], apData = factor(orderedPhyloData[, 3]))

  ## plot gene names
  orderedRowNamesMat <- data.frame(x = rep(0, length(orderedRowNames)),
                                   y = seq(0.5, length(orderedRowNames)-0.5, 1),
                                   fillName = orderedRowNames)
  

  geneNamesObj <- ggplot(orderedRowNamesMat, aes_string('x', 'y', label = 'fillName')) +
    geom_text(size = geneNameSize, colour = geneNameCol) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
            theme_phylo()


  ## plot phylogenetic matrix
  phyloObj <- ggplot(orderedPhyloData, aes_string('speNames', 'geneNames')) +
    geom_tile(aes_string(fill = 'apData')) +
      scale_fill_manual(name = 'status', labels = c('absent', 'present'), values = c(absentCol, presentCol)) +
        labs(x = NULL, y = NULL) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
            scale_x_continuous(expand = c(0, 0), breaks = NULL) +
              theme_phylo()

  ## dendrogram plot for genes
  ddata <- dendro_data.hclust(hcGene, type = 'rectangle')
  segData <- segment(ddata)
  segData[, c(1, 3)] <- segData[, c(1, 3)] - 0.5
  geneDendroObj <- ggplot(segData) +
    geom_segment(aes_string(x = 'x', y = 'y', xend = 'xend', yend = 'yend')) +
      labs(x = NULL, y = NULL) +
        scale_y_reverse(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(phyloData)), breaks = NULL) +
            coord_flip() +
              theme_phylo()

  ## gene color block
  orderedGeneColMat <- data.frame(x = rep(0, length(orderedGeneCol)),
                                  y = 1:length(orderedGeneCol),
                                  fillCol = factor(orderedGeneCol))

  geneBlockObj <- ggplot(orderedGeneColMat, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fillCol'), color = geneBetweenBlockCol) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedGeneColMat$fillCol)) +
              theme_phylo()

  ## species color block
  orderedSpeColMat <- data.frame(y = rep(0, length(orderedSpeCol)),
                                 x = 1:length(orderedSpeCol),
                                 fillCol = factor(orderedSpeCol))

  speBlockObj <- ggplot(orderedSpeColMat, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fillCol')) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedSpeColMat$fillCol)) +
              theme_phylo()
 
  
  ## plot empty block
  emptyData <- data.frame(x = 1, y = 1)
  empty <- ggplot(emptyData) +
    geom_point(aes_string('x', 'y'), colour='white') +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            theme_phylo()
  

  ## plotRes <- marrangeGrob(
  ##   list(empty, empty, empty, speBlockObj, geneDendroObj, geneNamesObj, geneBlockObj, phyloObj),
  ##   ncol = 4,
  ##   nrow = 2,
  ##   widths = widthsShinkage,
  ##   heights = heightsShinkage,
  ##   top = NULL)

  plotRes <- grid.arrange(
    empty,
    empty,
    empty,
    speBlockObj,
    geneDendroObj,
    geneNamesObj,
    geneBlockObj,
    phyloObj,
    ncol = 4,
    nrow = 2,
    widths = widthsShinkage,
    heights = heightsShinkage)

  return(plotRes)
} 



##' Phylo theme
##'
##' A totally blank gpplot2 theme
##' @title Blank theme for plots
##' @param ... ggplot2 theme() additional parameters
##' @return ggplot2 object
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 element_blank theme_bw %+replace% theme
##' @importFrom grid unit
##' 
theme_phylo <- function(...) {
  theme_bw() %+replace%
  theme(legend.position='none',
        title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, 'mm'),
        axis.ticks.margin = unit(0, 'mm'),
        axis.line = element_blank(),
        panel.margin = unit(0, 'mm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'line'),
        legend.margin = unit(0, 'mm'),
        ...)
}


##' Phylogenetic profiles of human F1Fo ATP synthase subunits.
##'
##' A list containing three objects:
##' "atpPhylo" is a numeric matrix representing the phylogenetic profiles of F1Fo ATP synthase subuints.
##' "specol" is the vector of colors corresponding to the species in "atpPhylo".
##' "genecol" is the vector of colors corresponding to the genes in "atpPhylo".
##' @docType data
##' @name fatp
##' @format A list
##' @references Unpublished data from Yulong Niu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' 
NULL
  
