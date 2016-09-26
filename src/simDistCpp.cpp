#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simDistCpp.h"

using namespace Rcpp;
using namespace arma;

//' @inheritParams SimCor
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double SimJaccard(arma::umat pairProfile) {

  arma::uvec f = pairProfile.col(0);
  arma::uvec t = pairProfile.col(1);

  arma::uvec binft = f + 2*t;

  arma::uvec inter = find(binft == 3);

  double A = inter.n_elem;

  double jac = A / (sum(f) + sum(t) - A);

  return jac;
}
