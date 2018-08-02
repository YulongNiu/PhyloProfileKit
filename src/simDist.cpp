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
double SimCor_(arma::vec f,
               arma::vec t) {
  mat corMat = cor(f, t);
  return corMat(0, 0);
}

//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double SimJaccard_(arma::vec f,
                   arma::vec t) {
  vec combVec = f + 2*t;
  double A = sum(combVec == 3);
  return  A / (sum(f) + sum(t) - A);
}


//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double DistHamming_(arma::vec f,
                    arma::vec t) {
  return sum(f != t);
}


//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double DistManhattan_(arma::vec f,
                      arma::vec t) {
  return sum(f != t);
}



//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double DistEuclidean_(arma::vec f,
                     arma::vec t) {
  return sqrt(sum(square(f - t)));
}


//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
// [[Rcpp::export]]
double DistMinkowski_(arma::vec f,
                      arma::vec t,
                      arma::uword p) {
  return pow(accu(pow(abs(f - t), p)), 1.0 / p);
}


