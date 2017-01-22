##' @include PP.R
NULL

##' The constructor of the \code{PP} class from the STRING database
##'
##' Retrieve the bit score matrix from a give organism. The default bit score threshold in STRING is 50, homologous with values smaller than 50 are not shown. For these homologous, the default value is set as 0.
##'
##' @title Constructor from the STRING database
##' @param speID Number. NCBI taxonomy identifier must be present in the STRING database. Use \code{get_STRING_species()} from the STRINGdb package to determine supported organisms.
##' @param idx A numeric vector or "all". 
##' @param STRINGnorm Whether you use the normalized or raw bit scores.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @param symbets Whether to keep only symmetrical best hits.
##' @param verbose Whether to print retrieve process.
##' @return A \code{PP} object.
##' @examples
##' ## human top 5 proteins
##' PPSTRING(9606, idx = 1:5, n = 2)
##'
##' \dontrun{
##' ## all human proteins
##' PPRCSTRING(9606, idx = "all", n = 4)
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @import STRINGdb
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom magrittr %>%
##' @importFrom methods new
##' @export
##' 
PPSTRING <- function(speID, idx = 1:10, symbets = FALSE, STRINGnorm = FALSE, n = 1, verbose = FALSE) {

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
  if (is.numeric(idx)) {
    idx <- idx[idx <= pIDsLen]
    if (length(idx) == 0) {
      return(new('PP'))
    } else {}
  } else {
    ## all proteins
    idx <- seq_len(pIDsLen)
  }

  proMat <- foreach(i = 1:length(idx), .combine = rbind) %dopar% {

    if (verbose) {
      cat(paste0('It is running ', i, ' in a total of ', length(idx), '.\n'))
    } else {}

    ## initvec
    eachInit <- rep(0, length(spIDs))

    eachM <- sdb$get_homologs_besthits(pIDs[idx[i]], symbets = symbets, bitscore_threshold = 0)

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
  rownames(proMat) <- pIDs[idx]
  colnames(proMat) <- spIDs

  ## stop multiple core
  stopImplicitCluster()

  proMat %>% PP %>% return

}
