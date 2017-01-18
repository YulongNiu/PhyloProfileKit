##' @include PP.R
NULL

##' Normalization methods
##'
##' Normalization of raw bit score profile by Singular Value Decomposition (SVD) or Normalized Phylogenetic Profile (NPP).
##'
##' The SVD is retrieved from the \href{https://bitbucket.org/andrea/svd-phy}{SVD-Phy package} with enhanced performance. The NPP is firstly described in the \href{http://www.nature.com/nature/journal/v493/n7434/full/nature11779.html}{paper}.
##'
##' @title Normalization phylogenetic profile
##' @inheritParams Norm
##' @return A \code{PP} object.
##' @examples
##' require("magrittr")
##' ppRawBit <- matrix(sample(0:200, 2000 * 1000, replace = TRUE),
##'                    ncol = 1000,
##'                    dimnames = list(paste0('protein', 1:2000),
##'                                    paste0('spe', 1:1000))) %>% PP
##'
##' ## SVD
##' Norm(ppRawBit, method = 'SVD', bitCutoff = 100, bitReset = 0, minConserve = 12, trimming = 0.8)
##'
##' ## NPP
##' Norm(ppRawBit, method = 'NPP', bitCutoff = 60, bitReset = 1, minConserve = 12)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \href{https://bitbucket.org/andrea/svd-phy}{SVD-Phy package}
##' @references \href{http://bioinformatics.oxfordjournals.org/content/suppl/2015/11/25/btv696.DC1/SVD-Phy-supplementary-material.docx}{SVD description}
##' @references \href{http://www.nature.com/nature/journal/v493/n7434/extref/nature11779-s1.pdf}{NPP description}
##' @rdname Norm-methods
##' @exportMethod Norm
##'
setMethod(f = 'Norm',
          signature = c(x = 'PP', method = 'character'),
          definition = function(x, method, ...) {

            d <- PPData(x)

            ## select norm function
            normd <- switch(method,
                            SVD = SVDNorm(d, ...),
                            NPP = NPPNorm(d, ...))

            PPData(x) <- normd
            return(x)
          })

