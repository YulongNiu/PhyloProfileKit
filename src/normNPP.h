#ifndef _NORMNPP_H_
#define _NORMNPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericMatrix NPPNorm(Rcpp::NumericMatrix rawBitM,
                            double bitCutoff,
                            double bitReset,
                            double minConserve);

#endif
