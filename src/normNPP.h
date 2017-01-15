#ifndef _NORMNPPCPP_H_
#define _NORMNPPCPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericMatrix NPPNorm(Rcpp::NumericMatrix rawBitM,
                            double bitCutoff,
                            double bitReset,
                            double minConserve);

#endif
