##' Plot correlation matrix of phylogenetic profiles.
##'
##' A combination plot of correlation matrix with genes and species clusters.
##' @title Plot correlation matrix
##' @param gradientCol The gradien colours for correlation matrix.
##' @param showCorVal Whether or not show the correlation values, and the default is set as TRUE.
##' @inheritParams PlotPhyloProfile
##' @importFrom grDevices colorRampPalette
##' @return A plot object.
##' @examples
##' data(fatp)
##' ATPCorPlot <- PlotPhyloCor(fatp$atpPhylo, geneCol = fatp$genecol)
##' \dontrun{
##' require(grDevices)
##' cairo_pdf('FATPCorplot.pdf')
##' ATPCorPlot <- PlotPhyloCor(fatp$atpPhylo, geneCol = fatp$genecol)
##' dev.off()
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_text geom_tile geom_segment geom_point scale_fill_manual labs scale_x_continuous scale_y_continuous scale_fill_gradientn scale_x_discrete scale_y_discrete scale_x_reverse scale_y_reverse aes_string coord_flip
##' @importFrom ggdendro dendro_data segment
##' @importFrom gridExtra grid.arrange
##' @importFrom reshape2 melt
##' @importFrom RColorBrewer brewer.pal
##' @rdname simplot
##' @export
##'
PlotPhyloCor <- function(phyloData,
                         gradientCol = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
                         geneNameSize = 3,
                         geneNameCol = 'grey55',
                         geneBetweenBlockCol = NA,
                         geneCol,
                         showCorVal = TRUE,
                         widthsShinkage = c(0.9, 0.7, 0.3, 7),
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
  corMelt <- corMelt[!(is.na(corMelt[, 3])),]
  corMelt <- data.frame(corMelt)
  colnames(corMelt) <- c('From', 'To', 'Cor')
  

  ## order 'geneCol'
  orderedGeneCol <- geneCol[match(rownames(orderedPhyloData), names(geneCol))]

  ## plot cor matrix
  simBreaks = seq(0, 1, 0.25)
  corMatObj <- ggplot(corMelt, aes_string('From', 'To', fill = 'Cor')) +
    geom_tile() +
      scale_fill_gradientn(colours = gradientCol, breaks = simBreaks, labels = format(simBreaks), limits = c(0, 1)) +
        scale_x_discrete(expand = c(0, 0), breaks = NULL) +
          scale_y_discrete(expand = c(0, 0), breaks = NULL) +
            labs(x = NULL, y = NULL) +
              theme_pp(legend.justification = c(1, 1),
                          legend.position = c(1, 1))

  if (showCorVal) {
    corMatObj <- corMatObj +
      geom_text(aes_string('From', 'To', label = 'Cor'), colour = '#073642', size = 4)} else {}

  ## ## plot legent
  ## ## Extract Legend
  ## corMatObjLegent <- ggplot(corMelt, aes_string('From', 'To', fill = 'Cor')) +
  ##   geom_tile() +
  ##     scale_fill_gradientn(colours = gradientCol)
  ## g_legend <- function(a.gplot){ 
  ##   tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  ##   leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box') 
  ##   legend <- tmp$grobs[[leg]] 
  ##   return(legend)}
  ## corMatLegent <- g_legend(corMatObjLegent)
  
  ## plot row gene names
  orderedRowNamesMat <- data.frame(x = rep(0, length(orderedRowNames)),
                                   y = seq(0.5, length(orderedRowNames)-0.5, 1),
                                   fillName = rev(orderedRowNames))
  

  geneRowNamesObj <- ggplot(orderedRowNamesMat, aes_string('x', 'y', label = 'fillName')) +
    geom_text(size = geneNameSize, colour = geneNameCol) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
            theme_pp(legend.position='none')


  ## plot col gene names
  orderedColNamesMat <- data.frame(y = rep(0, length(orderedRowNames)),
                                   x = seq(0.5, length(orderedRowNames)-0.5, 1),
                                   fillName = orderedRowNames)
  

  geneColNamesObj <- ggplot(orderedColNamesMat, aes_string('x', 'y', label = 'fillName')) +
    geom_text(size = geneNameSize, colour = geneNameCol, angle = 90) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
            theme_pp()

  
  ## dendrogram plot for genes
  ddata <- dendro_data(hcGene, type = 'rectangle')
  segData <- segment(ddata)
  segData[, c(1, 3)] <- segData[, c(1, 3)] - 0.5
  geneDendroObj <- ggplot(segData) +
    geom_segment(aes_string(x = 'x', y = 'y', xend = 'xend', yend = 'yend')) +
      labs(x = NULL, y = NULL) +
        scale_y_reverse(expand = c(0, 0), breaks = NULL) +
          scale_x_reverse(expand = c(0, 0), limits = c(nrow(phyloData), 0), breaks = NULL) +
            coord_flip() +
              theme_pp(legend.position='none')


  ## gene color block
  orderedGeneColMat <- data.frame(x = rep(0, length(orderedGeneCol)),
                                  y = 1:length(orderedGeneCol),
                                  fillCol = rev(factor(orderedGeneCol)))

  geneBlockObj <- ggplot(orderedGeneColMat, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fillCol'), color = geneBetweenBlockCol) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedGeneColMat$fillCol)) +
              theme_pp(legend.position='none')

  
  ## plot empty block
  emptyBlock <- pp_empty(colour = 'white')

  ## plotRes <- list(geneDendroObj = geneDendroObj,
  ##                 geneRowNamesObj = geneRowNamesObj,
  ##                 geneBlockObj = geneBlockObj,
  ##                 corMatObj = corMatObj,
  ##                 empty = empty,
  ##                 geneColNamesObj = geneColNamesObj,
  ##                 corMelt = corMelt)
  
  ## plotRes <- arrangeGrob(geneDendroOabj,
  ##                        geneRowNamesObj,
  ##                        geneBlockObj,
  ##                        corMatObj,
  ##                        empty,
  ##                        empty,
  ##                        empty,
  ##                        geneColNamesObj,
  ##                        ncol = 4,
  ##                        nrow = 2,
  ##                        widths = widthsShinkage,
  ##                        heights = heightsShinkage)

  plotRes <- grid.arrange(geneDendroObj,
                         geneRowNamesObj,
                         geneBlockObj,
                         corMatObj,
                         emptyBlock,
                         emptyBlock,
                         emptyBlock,
                         geneColNamesObj,
                         ncol = 4,
                         nrow = 2,
                         widths = widthsShinkage,
                         heights = heightsShinkage)
  return(plotRes)
  
}

##' @inheritParams PlotPhyloCor
##' @return correlation matrix
##' @examples
##' data(fatp)
##' corMat <- GetPhyloCorMat(fatp$atpPhylo)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom stats cor dist hclust
##' @rdname simplot
##' @export
GetPhyloCorMat <- function(phyloData) {
  ## cluster genes and species
  hcGene <- hclust(dist(phyloData), method = 'average')
  rowInd <- hcGene$order

  ## order 'phyloData'
  orderedPhyloData <- phyloData[rowInd, ]
  orderedRowNames <- rownames(orderedPhyloData)

  ## correlation matrix
  corMat <- round(cor(t(orderedPhyloData)), digits = 2)

  return(corMat)
}
