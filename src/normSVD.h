#ifndef _NORMSVDCPP_H_
#define _NORMSVDCPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericMatrix SVDPhy(Rcpp::NumericMatrix bitM,
                           double trimming,
                           double minConserve);

Rcpp::NumericMatrix SVDNor(Rcpp::NumericMatrix rawBitM,
                           double bitCutoff,
                           double bitReset,
                           double trimming,
                           double minConserve);

#endif
