#ifndef _SIMDISTCPP_H_
#define _SIMDISTCPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double SimJaccard(arma::umat pairProfile);

double SimMI(arma::umat pairProfile);

arma::uword DistHamming(arma::umat pairProfile);

double eachMI(double p1,
              double p2,
              double p3,
              double n);

#endif

