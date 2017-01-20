#ifndef _MICPP_H_
#define _MICPP_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double SimMIBin(arma::mat pairProfile);

double eachMI(double p1,
              double p2,
              double p3,
              double n);

double SimMIConti(arma::mat pairProfile,
                  arma::uword bin);

double Info(arma::uvec v,
            double n);

arma::uvec HistTwo(arma::vec x,
                   arma::vec y,
                   arma::uword bin);

arma::uvec FindInter(arma::vec x,
                     arma::vec interval);

arma::uword FindInterSingle(double value,
                            arma::vec interval);

arma::vec gInter (arma::vec x,
                  arma::uword bin);

arma::uvec CountRepeat (arma::uvec x);

#endif
