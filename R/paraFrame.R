## ##' Parallel framework for \code{PPIdx} and \code{PPTreeIdx}.
## ##'
## ##' \code{Batch()}: Batch process.
## ##'
## ##' @title Parallel framework
## ##' @param p 
## ##' @param ft 
## ##' @param n 
## ##' @param ... 
## ##' @return 
## ##' @examples 
## ##' @author Yulong Niu \email{niuylscu@@gmail.com}
## Batch <- function(p, ft, n = 1, FUN, ...) {

##   ## register multiple core
##   registerDoParallel(cores = n)

##   ppiNames <- rownames(ft)
##   ppiNum <- nrow(ft)

##   batchVec <- foreach(i = 1:ppiNum, .combine = c) %dopar% {
##     ## print(paste0('It is running ', i, ' in a total of ', ppiNum, '.'))
##     f <- t(p[ft[i, 1], ])
##     t <- t(p[ft[i, 2], ])
##     eachSD <- FUN(f, t, ...)

##     return(eachSD)
##   }

##   ## stop multiple core
##   stopImplicitCluster()

##   return(batchVec)

## }
