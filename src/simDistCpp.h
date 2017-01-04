#ifndef _SIMDISTCPP_H_
#define _SIMDISTCPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double SimCor(arma::mat pairProfile);

double SimJaccard(arma::mat pairProfile);

double SimMI(arma::mat pairProfile);

arma::uword DistHamming(arma::mat pairProfile);

double DistEuclidean(arma::mat pairProfile);

double eachMI(double p1,
              double p2,
              double p3,
              double n);

#endif

