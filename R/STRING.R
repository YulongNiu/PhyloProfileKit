##' @include PP.R
NULL

##' The constructor of the \code{PP} class from the STRING database
##'
##' Retrieve the bit score matrix from a give organism.
##'
##' @title Constructor from the STRING database
##' @param speID Number. NCBI taxonomy identifier must be present in the STRING database. Use \code{get_STRING_species()} from the STRINGdb package to determine supported organisms.
##' @param top Number or "all". Top number of proteins or the whole genome (set as "all")
##' @param STRINGnorm Whether you use the normalized or raw bit scores.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @param symbets Whether you keep only symmetrical best hits.
##' @return A \code{PP} object.
##' @examples
##' ## human top 5 proteins
##' PPSTRING(9606, top = 5, n = 2)
##'
##' \dontrun{
##' ## all human proteins
##' PPRCSTRING(9606, top = "all", n = 4)
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @import STRINGdb
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom magrittr %>%
##' @export
##' 
PPSTRING <- function(speID, top = 10, symbets = FALSE, STRINGnorm = FALSE, n = 1) {

  ## register multiple core
  registerDoParallel(cores = n)

  ## all species
  sp <- get_STRING_species(version = '10', species_name = NULL)
  spIDs <- sp[, 1]

  ## all proteins from a given organism
  sdb <- STRINGdb$new(version = '10', species = speID, score_threshold = 0)
  pMat <- sdb$get_proteins()
  pIDs <- pMat[, 1]
  pIDsLen <- length(pIDs)

  ## determine protein number
  if (is.numeric(top)) {
    if (top < pIDsLen) {
      tN <- top
    } else {
      tN <- pIDsLen
    }
  } else {
    ## all proteins
    tN <- pIDsLen
  }

  proMat <- foreach(i = 1:tN, .combine = rbind) %dopar% {

    ## initvec
    eachInit <- rep(0, length(spIDs))

    eachM <- sdb$get_homologs_besthits(pIDs[i], symbets = symbets, bitscore_threshold = 0)

    if (nrow(eachM) != 0) {

      if (STRINGnorm) {
        bs <- eachM[, 'best_hit_normscore']
      } else {
        bs <- eachM[, 'best_hit_bitscore']
      }

      eachInit[match(eachM[, 'species_id'], spIDs)] <- bs
    } else {}

    return(eachInit)

  }

  ## STRING IDs
  rownames(proMat) <- pIDs[1:tN]
  colnames(proMat) <- spIDs

  ## stop multiple core
  stopImplicitCluster()

  proMat %>% PP %>% return

}
