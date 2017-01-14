#ifndef _NORMSVDCPP_H_
#define _NORMSVDCPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat SVDArma(arma::mat bitM,
                  double trimming,
                  double minConserve);

#endif
