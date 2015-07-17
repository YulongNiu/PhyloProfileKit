##' Plot correlation matrix of phylogenetic profiles.
##'
##' A combination plot of correlation matrix with genes and species clusters.
##' @title Plot correlation matrix
##' @param gradientCol The gradien colours for correlation matrix.
##' @inheritParams
##' @return A plot object. 
##' @examples
##' data(atpPhyloExample)
##' load('../data/atpPhyloExample.RData')
##' ATPCorPlot <- plotPhyloCor(atpPhyloExample$atpPhylo, geneCol = atpPhyloExample$genecol)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_text geom_tile geom_segment scale_fill_manual labs scale_x_continuous scale_y_continuous scale_fill_gradientn scale_x_discrete scale_y_discrete theme aes element_blank coord_flip
##' @importFrom grid unit
##' @importFrom ggdendro dendro_data segment
##' @importFrom gridExtra grid.arrange
##' @importFrom reshape2 melt
##' @importFrom RColorBrewer brewer.pal
##' @export
##' 
plotPhyloCor <- function(phyloData,
                         gradientCol = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                         geneNameSize = 3,
                         geneNameCol = 'grey55',
                         geneBlockCol = NA,
                         geneCol,
                         widthShinkage = c(0.7, 0.7, 0.3, 7),
                         heightsShinkage = c(7, 0.7)) {
  ## require('grid')
  ## require('reshape2')
  ## require('ggplot2')
  ## require('RColorBrewer')
  ## require('ggdendro')
  ## require('gridExtra')

  ## cluster genes and species
  hcGene <- hclust(dist(phyloData), method = 'average')
  rowInd <- hcGene$order

  ## order 'phyloData'
  orderedPhyloData <- phyloData[rowInd, ]
  orderedRowNames <- rownames(orderedPhyloData)

  ## correlation matrix
  corMat <- round(cor(t(orderedPhyloData)), digits = 2)
  corMat[lower.tri(corMat)] <- NA
  corMat <- corMat[, ncol(corMat):1]
  corMelt <- melt(corMat)
  corMelt <- corMelt[-which(is.na(corMelt[, 3])),]
  corMelt <- data.frame(corMelt)
  colnames(corMelt) <- c('From', 'To', 'Cor')
  

  ## order 'geneCol'
  orderedGeneCol <- geneCol[match(rownames(orderedPhyloData), names(geneCol))]

  ## plot cor matrix
  corMatObj <- ggplot(corMelt, aes(From, To, fill = Cor)) +
    geom_tile() +
      geom_text(aes(From, To, label = Cor), color = '#073642', size = 4) +
  scale_fill_gradientn(colours = gradientCol) +
    scale_x_discrete(expand = c(0, 0), breaks = NULL) +
      scale_y_discrete(expand = c(0, 0), breaks = NULL) +
        labs(x = NULL, y = NULL) +
          theme(legend.justification = c(1, 1),
                legend.position = c(1, 1),
                title = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks.length = unit(0, "mm"),
                axis.ticks.margin = unit(0, "mm"),
                axis.line = element_blank(),
                panel.margin = unit(0, 'mm'),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.margin = unit(c(0, 0, 0, 0), 'line'),
                legend.margin = unit(0, 'mm'))

  ## ## plot legent
  ## ## Extract Legend
  ## corMatObjLegent <- ggplot(corMelt, aes(From, To, fill = Cor)) +
  ##   geom_tile() +
  ##     scale_fill_gradientn(colours = gradientCol)
    ## g_legend <- function(a.gplot){ 
  ##   tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  ##   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  ##   legend <- tmp$grobs[[leg]] 
  ##   return(legend)}
  ## corMatLegent <- g_legend(corMatObjLegent)
  
  ## plot row gene names
  orderedRowNamesMat <- data.frame(x = rep(0, length(orderedRowNames)),
                                   y = seq(0.5, length(orderedRowNames)-0.5, 1),
                                   fillName = rev(orderedRowNames))
  

  geneRowNamesObj <- ggplot(orderedRowNamesMat, aes(x, y, label = fillName)) +
    geom_text(size = geneNameSize, colour = geneNameCol) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
            theme(legend.position='none',
                  title = element_blank(),
                  axis.text = element_blank(),
                  axis.title = element_blank(),
                  axis.ticks.length = unit(0, "mm"),
                  axis.ticks.margin = unit(0, "mm"),
                  axis.line = element_blank(),
                  panel.margin = unit(0, 'mm'),
                  panel.grid = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.margin = unit(c(0, 0, 0, 0), 'line'),
                  legend.margin = unit(0, 'mm'))


  ## plot col gene names
  orderedColNamesMat <- data.frame(y = rep(0, length(orderedRowNames)),
                                   x = seq(0.5, length(orderedRowNames)-0.5, 1),
                                   fillName = orderedRowNames)
  

  geneColNamesObj <- ggplot(orderedColNamesMat, aes(x, y, label = fillName)) +
    geom_text(size = geneNameSize, colour = geneNameCol, angle = 90) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
            theme(legend.position='none',
                  title = element_blank(),
                  axis.text = element_blank(),
                  axis.title = element_blank(),
                  axis.ticks.length = unit(0, "mm"),
                  axis.ticks.margin = unit(0, "mm"),
                  axis.line = element_blank(),
                  panel.margin = unit(0, 'mm'),
                  panel.grid = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.margin = unit(c(0, 0, 0, 0), 'line'),
                  legend.margin = unit(0, 'mm'))
  
  ## dendrogram plot for genes
  ddata <- dendro_data(hcGene, type = "rectangle")
  segData <- segment(ddata)
  segData[, c(1, 3)] <- segData[, c(1, 3)] - 0.5
  geneDendroObj <- ggplot(segData) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      labs(x = NULL, y = NULL) +
        scale_y_reverse(expand = c(0, 0), breaks = NULL) +
          scale_x_reverse(expand = c(0, 0), limits = c(nrow(phyloData), 0), breaks = NULL) +
            coord_flip() +
              theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))

  ## gene color block
  orderedGeneColMat <- data.frame(x = rep(0, length(orderedGeneCol)),
                                  y = 1:length(orderedGeneCol),
                                  fillCol = rev(factor(orderedGeneCol)))

  geneBlockObj <- ggplot(orderedGeneColMat, aes(x, y)) +
    geom_tile(aes(fill = fillCol), color = geneBlockCol) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedGeneColMat$fillCol)) +
              theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))
  
   ## plot empty block
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    labs(x = NULL, y = NULL) +
      scale_y_continuous(expand = c(0, 0), breaks = NULL) +
        scale_x_continuous(expand = c(0, 0), breaks = NULL) +
          theme(legend.position='none',
                title = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks.length = unit(0, "mm"),
                axis.ticks.margin = unit(0, "mm"),
                axis.line = element_blank(),
                panel.margin = unit(0, 'mm'),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.margin = unit(c(0, 0, 0, 0), 'line'),
                legend.margin = unit(0, 'mm'))
  
  plotRes <- grid.arrange(geneDendroObj, geneRowNamesObj, geneBlockObj, corMatObj, empty, empty, empty, geneColNamesObj, ncol = 4, nrow = 2, widths = widthShinkage, heights = heightsShinkage)
  
  return(plotRes)
  
}
