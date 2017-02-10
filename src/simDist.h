#ifndef _SIMDIST_H_
#define _SIMDIST_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double SimCor(arma::vec f,
              arma::vec t);

double SimJaccard(arma::vec f,
                  arma::vec t);

double DistHamming(arma::vec f,
                   arma::vec t);

double DistEuclidean(arma::vec f,
                     arma::vec t);

#endif

