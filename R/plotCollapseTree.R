## htMat <- sceTreeL@.Data[c(1, 50), ]
## htMatChar <- apply(htMat, 0:1, as.character)
## rownames(htMatChar) <- colnames(htMat)

## p <- ggtree(sceTree) + geom_tiplab()

## gheatmap(p, htMatChar, width = 0.5) +
##     scale_fill_manual(breaks = c('0', '1'), values=c('gray', 'steelblue'))
