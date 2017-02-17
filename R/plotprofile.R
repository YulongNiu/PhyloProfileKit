##' The core plot of phylogenetic profile
##'
##' Plot of phylogenetic profile with protein names, a protein colour bar (optional), and a species colour bar (optional).
##'
##' @title Core plot of phylogenetic profile
##' @param p A integer matrix (binning profile).
##' @param proSize A numeric value. Size of protein names.
##' @param proGroup A factor indicating protein groups.
##' @param proGroupCol A named character vector indicating protein group colour.
##' @param speGroup A factor indicating species groups.
##' @param speGroupCol A named character vector indicating species group colour.
##' @return A \code{gmat} object.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 scale_fill_manual coord_flip scale_fill_gradientn
##' @keywords internal
##' 
ProfileCore <- function(p,
                        proSize = 3,
                        proGroup = NA,
                        proGroupCol = NA,
                        speGroup = NA,
                        speGroupCol = NA){

  ## define color
  binColor <- c('0' = '#eff3ff', '1' = '#2171b5')
  contiColor <- colorRampPalette(c('#eff3ff', '#bdd7e7', '#6baed6', '#2171b5'))(100)

  ## center
  cObj <- pp_profile(p)
  if (isBinMat_internal(p)) {
    cObj <- cObj + scale_fill_manual(values = binColor)
  } else {
    cObj <- cObj + scale_fill_gradientn(colours = contiColor, breaks = seq(min(p), max(p), length.out = 5))
  }

  ## left
  if (!is.na(proGroup)) {
    lObj <- pp_tile(proNames) +
      scale_fill_manual(values = proGroupCol)
  } else {
    lObj <- list()
  }
  lObj <- lObj %@+% pp_text(rownames(p), size = proSize)

  ## top
  if (!is.na(speGroup)) {
    tObj <- pp_tile(speGroup) +
      scale_fill_manual(values = speGroupCol) +
      coord_flip()
  } else {
    tObj <- list()
  }

  pObj <- ascore(cObj) %@<% lObj %@^% tObj

  return(pObj)
}
