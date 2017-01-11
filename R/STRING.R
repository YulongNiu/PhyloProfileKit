##' @include PP.R
NULL

##' The constructor of the \code{PP} class from the STRING database
##'
##' Retrieve the bit score matrix from a give organism.
##'
##' @title Constructor from the STRING database
##' @param speID Number. NCBI taxonomy identifier must be present in the STRING database.
##' @param STRINGnorm Whether use the normalized or raw bit scores.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @inheritParams STRINGdb::get_homologs_besthits
##' @return A \code{PP} object.
##' @examples
##' ## human
##' PPSTRING(9606, n = 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom STRINGdb get_STRING_species get_proteins get_proteins
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

  proMat <- foreach(i = 1:10, .combine = rbind) %dopar% {

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

  ## CAUTIONS!!!!
  rownames(proMat) <- pIDs[1:10]
  colnames(proMat) <- spIDs

  ## stop multiple core
  stopImplicitCluster()

  proMat %>% PP %>% return

}
