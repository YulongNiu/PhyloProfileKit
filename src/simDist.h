#ifndef _SIMDIST_H_
#define _SIMDIST_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double SimCor(arma::mat pairProfile);

double SimJaccard(arma::mat pairProfile);

arma::uword DistHamming(arma::mat pairProfile);

double DistEuclidean(arma::mat pairProfile);

#endif

