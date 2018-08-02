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

double DistManhattan_(arma::vec f,
                      arma::vec t);

double DistEuclidean_(arma::vec f,
                      arma::vec t);

double DistMinkowski_(arma::vec f,
                      arma::vec t,
                      arma::uword p);

#endif

