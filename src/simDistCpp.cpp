#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simDistCpp.h"

using namespace Rcpp;
using namespace arma;

//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double SimJaccard(arma::umat pairProfile) {

  uvec f = pairProfile.col(0);
  uvec t = pairProfile.col(1);
  uvec combVec = f + 2*t;

  uvec inter = find(combVec == 3);

  double A = inter.n_elem;

  double jac = A / (sum(f) + sum(t) - A);

  return jac;
}


//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double SimMI(arma::umat pairProfile) {

  uvec f = pairProfile.col(0);
  uvec t = pairProfile.col(1);
  uvec combVec = f + 2*t;

  double N = pairProfile.n_rows;
  uvec Avec = find(combVec == 3);
  double A = Avec.n_elem;
  uvec Bvec = find(combVec == 1);
  double B = Bvec.n_elem;
  uvec Cvec = find(combVec == 2);
  double C = Cvec.n_elem;
  double D = N - A - B - C;

  double I = eachMI(A, B, C, N) + eachMI(B, A, D, N) + eachMI(C, A, D, N) + eachMI(D, C, B, N);

  return I;

}


//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
arma::uword DistHamming(arma::umat pairProfile) {

  uvec f = pairProfile.col(0);
  uvec t = pairProfile.col(1);

  uvec neq = find(f != t);

  uword ham = neq.n_elem;

  return ham;

}

//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double DistEuclidean(arma::umat pairProfile) {

  uvec f = pairProfile.col(0);
  uvec t = pairProfile.col(1);

  uvec neq = f - t;

  double eu = sqrt(sum(square(neq)));

  return eu;

}


//' @keywords internal
// [[Rcpp::export]]
double eachMI(double p1,
              double p2,
              double p3,
              double n) {
  if (p1 == 0) {
    return 0.0;
  } else {
    return p1 * log(n * p1 / ((p1 + p2) * (p1 + p3))) / n;
  }
}

