##' Extract segment data
##'
##' @title Extract segments
##' @param d Data structure from tree from the ggtree package.
##' @return A data frame.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
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


##' Deal with \code{phylo} class objects for plot
##'
##' \code{Phylo2Mat()}: x and y locations.
##'
##' \code{Phylo2MatX()}: x locations.
##'
##' \code{Phylo2MatY()}: y locations.
##'
##' \code{InitM()}: edge and edge length.
##'
##' \code{OrderedTip()}: Ordered tips indices in the tree plot.
##'
##' Roles of \code{phylo} class:
##'
##' \itemize{
##'   \item \code{tree$edge}: a two-column matrix, 1st --> 2nd is from root to tips.
##'   \item \code{edge.length}: same order with \code{tree$edge}.
##'   \item \code{tip.label}: The order is from \code{1:Ntip(tree)}
##'   \item In plot, the tips appear in the same order in \code{tree$edge}.
##' }
##'
##' @title tree utilities
##' @param x A \code{phylo} class
##' @return
##'
##' \code{Phylo2Mat()}: A data.frame indicating x and y locations.
##'
##' \code{Phylo2MatX()}: A data.frame indicating x locations.
##'
##' \code{Phylo2MatY()}: A data.frame indicating y locations.
##'
##' \code{InitM()}: A numeric matrix indicating edge and edge length.
##'
##' \code{OrderedTip()}: A numeric vector indicating the ordered tips indices in the tree plot.
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname tree-uti
##' @keywords internal
##' 
Phylo2Mat <- function(x) {
  locX <- Phylo2MatX(x)
  locY <- Phylo2MatY(x)

  m <- merge(locX, locY, by.x = 'node', by.y = 'node')

  return(m)
}


##' @inheritParams Phylo2Mat
##' @importFrom magrittr  %<>%
##' @rdname tree-uti
##' @keywords internal
##' 
Phylo2MatY <- function(x) {

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
  mutM %<>% as.data.frame

  return(mutM)
}

##' @inheritParams Phylo2Mat
##' @importFrom ape Ntip
##' @importFrom magrittr  %<>% %>%
##' @rdname tree-uti
##' @keywords internal
##' 
InitM <- function(x) {

  edgeMat <- x$edge
  tipNum <- Ntip(x)

  edgeMat %<>% rbind(rep(tipNum + 1, 2))

  initX <- rep(0, nrow(edgeMat))
  tipIdx <- edgeMat[, 2] <= tipNum
  initX[tipIdx] <- 1:tipNum


  edgeMat %<>% cbind(initX)

  return(edgeMat)
}


##' @inheritParams Phylo2Mat
##' @importFrom ape Ntip
##' @rdname tree-uti
##' @keywords internal
##' 
OrderedTip <- function(x) {
  tipNum <- Ntip(x)

  m <- InitM(x)[, 2]
  otIdx <- m[m <= tipNum]

  return(otIdx)
}


##' @inheritParams Phylo2Mat
##' @importFrom magrittr  %<>%
##' @rdname tree-uti
##' @keywords internal
##' 
Phylo2MatX <- function(x) {

  edgeMat <- x$edge
  edge <- x$edge.length

  start <- root <- Ntip(x) + 1
  startLen <- 0

  while(TRUE) {
    starttmp <- NULL
    startLentmp <- NULL

    for (i in seq_along(start)){
      eachEndIdx <- edgeMat[, 1] == start[i]
      if (sum(eachEndIdx) > 1) {
        edge[eachEndIdx] <- edge[eachEndIdx] + startLen[i]
        starttmp %<>% c(edgeMat[eachEndIdx, 2])
        startLentmp %<>% c(edge[eachEndIdx])
      } else {}
    }

    if (is.null(starttmp)) {
      break
    } else {
      start <- starttmp
      startLen <- startLentmp
    }
  }

  locX <- rbind(cbind(edgeMat, edge),
                c(root, root, 0))[, 2:3]

  colnames(locX) <- c('node', 'x')
  locX %<>% as.data.frame

  return(locX)
}
