##' Plot trees
##'
##' Plot phylogenetic trees and dendrograms
##'
##' @title pp plot trees and dendrograms
##' @param x A \code{phylo} or a \code{hclust} data.
##' @param ... Parameters passed to \code{geom_segment()} in the ggplot2 package.
##' @return A \code{gg} class object
##' @examples
##' require('ggplot2')
##' require('ape')
##'
##' treePath <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfile")
##' sceTree <- read.nexus(treePath)
##'
##' pp_tree(sceTree)
##' pp_tree(sceTree) + coord_flip() + scale_x_reverse()
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_tile labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @importFrom ape as.phylo Ntip
##' @seealso \code{\link[ggplot2]{geom_segment}}
##' @export
##' 
pp_tree <- function(x, ...) {

  if (inherits(x, 'hclust')) {
    x <- as.phylo(x)
} else {}

  segData <- ExtractSeg(Phylo2Mat(x))

  segData[, c(3, 4)] <- segData[, c(3, 4)] - 0.5

  tObj <- ggplot(segData) +
    geom_segment(aes_string(x = 'x', y = 'y', xend = 'xend', yend = 'yend'), ...) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, Ntip(x)), breaks = NULL) +
    theme_pp(legend.position='none')

  return(tObj)
}

##' Extract segment data
##'
##' @title Extract segments
##' @param d Data structure from tree from the ggtree package.
##' @return A data frame.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##' 
ExtractSeg <- function(d) {
  p <- d$parent

  seg1 <- data.frame(x = d$x[p],
                     xend = d$x,
                     y = d$y,
                     yend = d$y)

  seg2 <- data.frame(x = d$x[p],
                     xend = d$x[p],
                     y = d$y[p],
                     yend = d$y)

  seg <- rbind(seg1, seg2)

  return(seg)
}

require('PhyloProfile')
require('ggplot2')
require('ape')
require('ggtree')
require('magrittr')
treePath <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfile")
sceTree <- read.nexus(treePath)

pp_tree(sceTree) %@+% (ggtree(sceTree) + geom_tiplab(size=3, color="purple"))  %@+% (pp_tree(sceTree) + coord_flip() + scale_x_reverse()) %@+% pp_text(sceTree$tip.label)

pp_tree(sceTree) %@+% pp_text(sceTree$tip.label) %@+% (ggtree(sceTree) + geom_tiplab(size=3, color="purple"))

ggtree(sceTree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()


tmp1 <- rtree(10)
plot(tmp1, edge.width = 2)
nodelabels()
tiplabels()

Phylo2Mat <- function(x) {

}

Phylo2MatY <- function(x) {

  initM <- InitM(x)[, c(1, 2, 4)]

  locYList <- lapply(nodepath(x), EachBranch, initM)

  locY <- do.call(rbind, locYList)
  locY <- locY[!duplicated(locY), ]

  locY <- rbind(locY, c(Ntip(x) + 1, 0))
  colnames(locY) <- c('node', 'y')

  return(locY)
}

Phylo2MatX <- function(x) {

  initM <- InitM(x)[, 1:3]

  ## initiation
  mutM <- matrix(nrow = 0, ncol = ncol(initM))

  while(TRUE) {

    if (nrow(initM) == 0) {
      break;
    } else {}

    pMat <- initM[initM[, 3] != 0, ]
    pUni <- unique(pMat[, 1])

    for (i in pUni) {
      eachIdx <- pMat[, 1] == i
      if (sum(eachIdx) > 1) {
        loc <- mean(pMat[eachIdx, 3])
        initM[initM[, 2] == i, 3] <- loc

        moveIdx <- initM[, 1] == i
        mutM <- rbind(mutM, initM[moveIdx, ])

        initM <- initM[!moveIdx, ]
      } else {}
    }
  }

  colnames(mutM) <- c('parent', 'node', 'y')

  return(mutM)
}

InitM <- function(x) {

  edgeMat <- x$edge
  tipNum <- Ntip(x)
  edge <- x$edge.length

  edgeMat <- rbind(edgeMat, rep(tipNum + 1, 2))

  initX <- rep(0, nrow(edgeMat))
  tipIdx <- edgeMat[, 2] <= tipNum
  initX[tipIdx] <- edgeMat[tipIdx, 2]


  edgeMat %<>% cbind(initX) %>% cbind(c(edge, 0))

  return(edgeMat)
}


EachBranch <- function(x, m) {

  edgeNum <- length(x) - 1

  initY <- numeric(edgeNum)

  for (i in 1:edgeNum) {
    eachEdge <- c(x[i], x[i+1])
    eachIdx <- which((m[, 1] == eachEdge[1]) &
                     (m[, 2] == eachEdge[2]))
    initY[i] <- m[eachIdx, 3]
  }

  initY %<>% cumsum %>% cbind(x[2:length(x)], .)

  return(initY)
}
