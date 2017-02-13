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


##' Deal with \code{phylo} class objects for plot
##'
##' \code{Phylo2Mat()}: x and y locations.
##'
##' \code{Phylo2MatX()}: x locations.
##'
##' \code{Phylo2MatY()}: y locations.
##'
##' \code{EachBranch()}: y locations for each path.
##'
##' \code{InitM()}: edge and edge length.
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
##' \code{InitM()}: A numeric matrix indicating edge and edge length
##'
##' \code{EachBranch()}: A numeric matrix indicating y locations for each path.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
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
##' @importFrom ape Ntip nodepath
##' @rdname tree-uti
##' @keywords internal
##' 
Phylo2MatX <- function(x) {

  initM <- InitM(x)[, c(1, 2, 4)]

  locXList <- lapply(nodepath(x), EachBranch, initM)

  locX <- do.call(rbind, locXList)
  locX <- locX[!duplicated(locX), ]

  locX %<>% rbind(c(Ntip(x) + 1, 0))
  colnames(locX) <- c('node', 'x')

  locX %<>% as.data.frame

  return(locX)
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
  edge <- x$edge.length

  edgeMat %<>% rbind(rep(tipNum + 1, 2))

  initX <- rep(0, nrow(edgeMat))
  tipIdx <- edgeMat[, 2] <= tipNum
  initX[tipIdx] <- edgeMat[tipIdx, 2]


  edgeMat %<>% cbind(initX) %>% cbind(c(edge, 0))

  return(edgeMat)
}


##' @param x A numeric vector for a single path.
##' @param m A numeric matrix indicating edge length
##' @importFrom magrittr  %<>% %>%
##' @rdname tree-uti
##' @keywords internal
##' 
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
