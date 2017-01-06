##' SVD normalization
##'
##' Algorithm:
##'
##' Step1: rawBitM < bitCutoff to bitReset;
##'
##' Step2: In each row (species), x/max(x);
##'
##' Step3: L^2 SVD normalization.
##'
##' The core SVD normalization is retrieved from the SVD-Phy package with performance modification.
##' @title Singular value decomposition normalization of bit score matrix
##' @param rawBitM Raw bit score matrix.
##' @param bitCutoff Minimum value of the bit score.
##' @param bitReset Reset the bit score for ones lower than the `bitCutoff`.
##' @return
##'
##' SVDNor(): SVD normalized bit score matrix.
##'
##' SVDPhy(): A L^2 normalized unitary matrix.
##'
##' @examples
##' data(fatp)
##' svdM <- SVDNor(fatp$atpPhyloBit)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @inheritParams SVDPhy
##' @references \url{http://bioinformatics.oxfordjournals.org/content/suppl/2015/11/25/btv696.DC1/SVD-Phy-supplementary-material.docx}
##' @references \url{https://bitbucket.org/andrea/svd-phy}
##' @seealso \code{\link{NPPNor}}
##' @rdname SVD
##' @export
SVDNor <- function(rawBitM, bitCutoff = 60, bitReset = 1, trimming = 0.3, minConserve = -0.1) {

  ## step1: rawBitM < bitCutoff to bitReset
  norProfile <- apply(rawBitM, 1:2, function(x){
    x <- ifelse(x < bitCutoff, bitReset, x)
    return(x)
  })

  ## step2: In each row (species), x/max(x)
  norProfile <- apply(norProfile, 1, function(x) {
    x <- x/max(x)
    return(x)
  })
  norProfile <- t(norProfile)

  ## step3: L^2 SVD normalization.
  norProfile <- SVDPhy(norProfile, trimming = trimming, minConserve = minConserve)

  return(norProfile)

}




##' @param bitM Bit score matrix, for example the BLASTP or STRING bit scores. It is a named numeric matrix, columns are species and rows are genes.
##' @param trimming A percentages top unitary matrix.
##' @param minConserve Minimum number of homologous in each species. The species with homologous less than thsi value are discarded.
##' @rdname SVD
##' @export
SVDPhy <- function(bitM, trimming, minConserve){

  ## Lapack SVD
  s <- La.svd(bitM)
  svdm <- s$u

  ## filter genes
  filteredIdx <- apply(bitM, 1, function(x){
    counter <- sum(x > 0)
    return(counter > minConserve)
  })
  svdm <- svdm[filteredIdx, ]

  ## trim species
  svdmTrimmed <- svdm[, 1:round(trimming * ncol(svdm))]

  ## L^2 normalization
  resultM <- t(apply(svdmTrimmed, 1, function(x){
    return(x/sqrt(sum(x^2)))
  }))

  rownames(resultM) <- rownames(bitM)

  return(resultM)
}
