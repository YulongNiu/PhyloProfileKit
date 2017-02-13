##' Plot trees
##'
##' Plot phylogenetic trees and dendrograms
##'
##' @title pp plot trees and dendrograms
##' @param x A \code{phylo} or a \code{hclust} data.
##' @param shift A numeric value indicating shit scale.
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
##' pp_tree(sceTree) %@+% pp_text(sceTree$tip.label)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_tile labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @importFrom ape as.phylo Ntip
##' @seealso \code{\link[ggplot2]{geom_segment}}
##' @export
##' 
pp_tree <- function(x, shift = -0.5, ...) {

  if (inherits(x, 'hclust')) {
    x <- as.phylo(x)
} else {}

  segData <- ExtractSeg(Phylo2Mat(x))

  segData[, c(3, 4)] <- segData[, c(3, 4)] + shift

  tObj <- ggplot(segData) +
    geom_segment(aes_string(x = 'x', y = 'y', xend = 'xend', yend = 'yend'), ...) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, Ntip(x)), breaks = NULL) +
    theme_pp(legend.position='none')

  return(tObj)
}

