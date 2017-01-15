##' NPP normalization
##'
##' Algorithm:
##'
##' Step1: rawBitM < hitCutoff to hitReset;
##'
##' Step2: 
##'
##' Step3: in each row (species), log2(x/max(x));
##'
##' Step4: z-score for each column.
##' 
##' @title z-score normalization of phylogenetic profile
##' @return
##'
##' NPPNor(): NPP normalized bit score matrix.
##'
##' @examples
##' data(fatp)
##' nppM <- NPPNorR(fatp$atpPhyloBit)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @inheritParams SVDNorR
##' @references \url{http://www.nature.com/nature/journal/v493/n7434/extref/nature11779-s1.pdf}
##' @seealso \code{\link{SVDNorR}}
##' @export
##' 
NPPNorR <- function(rawBitM, bitCutoff = 50, bitReset = 1) {

  # Step1: rawBitM < hitCutoff to hitReset;
  norProfile <- apply(rawBitM, 1:2, function(x){
    x <- ifelse(x < bitCutoff, bitReset, x)
    return(x)
  })

  ## step2: in each row (species), log2(x/max(x))
  norProfile <- apply(norProfile, 1, function(x) {
    x <- log2(x/max(x))
    return(x)
  })
  norProfile <- t(norProfile)

  ## step3: z-score for each column
  norProfile <- apply(norProfile, 2, function(x) {
    x <- scale(x)
    return(x)
  })

  rownames(norProfile) <- rownames(rawBitM)

  return(norProfile)
}
