#ifndef _NORMSVDCPP_H_
#define _NORMSVDCPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericMatrix SVDNor(Rcpp::NumericMatrix rawBitM,
                           double bitCutoff,
                           double bitReset,
                           double minConserve,
                           double trimming);

Rcpp::NumericMatrix SVDPhy(Rcpp::NumericMatrix bitM,
                           double bitReset,
                           double minConserve,
                           double trimming);

#endif
