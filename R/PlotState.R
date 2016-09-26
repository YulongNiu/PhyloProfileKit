NodeState <- function(tree, nodePath, pr) {
  
  fromEdge <- c(PhyloProTree:::InferEdge(tree$edge, nodePath, pr))
  names(fromEdge) <- c(tree$edge)
  fromEdge <- fromEdge[!duplicated(c(tree$edge))]

  fromEdge <- fromEdge[order(as.numeric(names(fromEdge)))]

  return(fromEdge)
}

CombineNodeState <- function(tree, nodePath, pr1, pr2) {
  
  state1 <- NodeState(tree, nodePath, pr1)
  state2 <- NodeState(tree, nodePath, pr2)

  state <- cbind(state1, state2)

  return(state)
}

PlotState <- function(tree, cState) {
  cState <- apply(cState, 1, paste, collapse = ',')

  state <- lapply(cState, function(x) {
    eachState <- data.frame(label = x)
    ggplot(eachState, aes(x = 0, y = 0, label = label)) + geom_text() + theme_inset()
  })
  names(state) <- 1:length(state)
  
  p <- ggtree(tree) + geom_tiplab(color="purple", hjust= -0.2)
  inset(p, state)
}

