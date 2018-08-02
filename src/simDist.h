#ifndef _SIMDIST_H_
#define _SIMDIST_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double SimCor_(arma::vec f,
               arma::vec t);

double SimJaccard_(arma::vec f,
                   arma::vec t);

double DistHamming_(arma::vec f,
                    arma::vec t);

double DistEuclidean_(arma::vec f,
                      arma::vec t);

#endif

