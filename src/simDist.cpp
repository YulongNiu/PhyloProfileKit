#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simDist.h"

using namespace Rcpp;
using namespace arma;

//' Similarity or distance of paired phylogenetic profile
//'
//' \code{SimCor()}: Person's correlation coefficient.
//'
//' \code{SimJaccard()}: Jaccard similarity.
//'
//' \code{SimMIBin()}: Mutual information for binning data.
//'
//' \code{SimMIConti()}: Mutual information for continuous data.
//'
//' \code{DistHamming()}: Hamming distance.
//'
//' \code{DistEuclidean()}: Euclidean distance.
//'
//' @title Similarity and distance
//' @param f Numeric vector indicating a gene profile.
//' @param t Numeric vector indicating a gene profile.
//' @return A numeric value.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double SimCor(arma::vec f,
              arma::vec t) {

  mat corMat = cor(f, t);
  double corVal = corMat(0, 0);

  return corVal;
}

//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double SimJaccard(arma::vec f,
                  arma::vec t) {

  vec combVec = f + 2*t;

  double A = sum(combVec == 3);

  double jac = A / (sum(f) + sum(t) - A);

  return jac;
}

//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double DistHamming(arma::vec f,
                   arma::vec t) {

  // Hamming distance
  // uword ham = sum(f != t);

  // Manhattan distance
  double ham = sum(abs(f - t));

  return ham;

}

//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double DistEuclidean(arma::vec f,
                     arma::vec t) {

  vec neq = f - t;

  double eu = sqrt(sum(square(neq)));

  return eu;
}

