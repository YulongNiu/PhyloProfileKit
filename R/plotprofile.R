##' @include AllClasses.R AllGenerics.R tree_utilities.R operators.R
NULL

##' Plot phylogenetic profile with additional data
##'
##' \itemize{
##'   \item x is a \code{PP} or \code{PPIdx} object: \code{method} works for both dimensions (proteins and species).
##'   \item x is a \code{PPTreeIdx} object: \code{method} works only for rows (proteins).
##' }
##'
##' @inheritParams plotprofile
##' @title Wrapped profiles plot function.
##' @return A \code{gtable} object.
##' @examples
##' data(fatp)
##'
##' p <- plotprofile(PP(fatp$atpPhylo),
##'                  method = 'euclidean',
##'                  proGroup = fatp$complex,
##'                  proGroupCol = fatp$complexCol,
##'                  speGroup = fatp$domain,
##'                  speGroupCol = fatp$domainCol)
##' 
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname plotprofile-methods
##' @importFrom magrittr %>%
##' @importFrom ape as.phylo
##' @exportMethod plotprofile
##'
setMethod(f = 'plotprofile',
          signature = c(x = 'PP'),
          definition = function(x, method, ...) {
            p <- x@.Data

            if (is.na(method)) {
              pObj <- ProfileCore(p, ...)
            } else {
              hcPro <- p %>% dist(method = method) %>% hclust(method = 'average')
              hcSpe <- p %>% t %>% dist(method = method) %>% hclust(method = 'average')

              ## the rownames(p) and hcPro$tip.label are the same
              hcOrder <- hcPro %>% as.phylo %>% OrderedTip
              pObj <- ProfileCore(p, rowIdx = hcOrder, colIdx = hcSpe$order, ...) %@<% pp_tree(hcPro)
            }

            ## pre-define widths and heights
            preWidths <- c(0.9, 0.7, 0.3, 7)
            preHeights <- c(0.3, 7)
            w <- rev(rev(preWidths)[seq_along(c(pObj@left, pObj@core))])
            h <- rev(rev(preWidths)[seq_along(c(pObj@top, pObj@core))])

            plotObj <- plot(pObj,
                            widths = w,
                            heights = h)

            return(plotObj)
          })



##' The core plot of phylogenetic profile
##'
##' Plot of phylogenetic profile with protein names, a protein colour bar (optional), and a species colour bar (optional).
##'
##' @title Core plot of phylogenetic profile
##' @param p A integer matrix (binning profile).
##' @param rowIdx A numeric vector indication row indices.
##' @param colIdx A numeric vector indication column indices.
##' @param proSize A numeric value. Size of protein names.
##' @param proGroup A factor indicating protein groups.
##' @param proGroupCol A named character vector indicating protein group colour.
##' @param speGroup A factor indicating species groups.
##' @param speGroupCol A named character vector indicating species group colour.
##' @return A \code{gmat} object.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom ggplot2 scale_fill_manual coord_flip scale_fill_gradientn
##' @keywords internal
##' 
ProfileCore <- function(p,
                        rowIdx = 1:nrow(p),
                        colIdx = 1:ncol(p),
                        proSize = 3,
                        proGroup = NA,
                        proGroupCol = NA,
                        speGroup = NA,
                        speGroupCol = NA){

  ## define color
  ## blues from colorbrewer2.org
  blue4 <- c('#eff3ff', '#bdd7e7', '#6baed6', '#2171b5')
  binColor <- c('0' = blue4[1], '1' = blue4[4])
  contiColor <- colorRampPalette(blue4)(100)

  ## rearrange params
  p <- p[rowIdx, colIdx]
  proGroup <- proGroup[rowIdx]
  speGroup <- speGroup[colIdx]

  ## center
  cObj <- pp_profile(p)
  if (isBinMat_internal(p)) {
    cObj <- cObj + scale_fill_manual(values = binColor)
  } else {
    cObj <- cObj + scale_fill_gradientn(colours = contiColor, breaks = seq(min(p), max(p), length.out = 5))
  }

  ## left
  if (!all(is.na(proGroup))) {
    lObj <- pp_tile(proGroup) +
      scale_fill_manual(values = proGroupCol)
  } else {
    lObj <- list()
  }
  lObj <- lObj %@+% pp_text(rownames(p), size = proSize)

  ## top
  if (!all(is.na(speGroup))) {
    tObj <- pp_tile(speGroup) +
      scale_fill_manual(values = speGroupCol) +
      coord_flip()
  } else {
    tObj <- list()
  }

  pObj <- ascore(cObj) %@<% lObj %@^% tObj

  return(pObj)
}
