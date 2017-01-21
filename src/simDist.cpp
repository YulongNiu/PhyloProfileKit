#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simDist.h"

using namespace Rcpp;
using namespace arma;

//' Similarity or distance of paired phylogenetic profile
//'
//' \code{SimCor()}: Person's correlation coefficient.
//' \code{SimJaccard()}: Jaccard similarity.
//' \code{SimMIBin()}: Mutual information for binning data.
//' \code{SimMIConti()}: Mutual information for continuous data.
//' \code{DistHamming()}: Hamming distance.
//' \code{DistEuclidean()}: Euclidean distance.
//'
//' @title similarity and distance
//' @param f Numeric vector indicating a gene profile.
//' @param t Numeric vector indicating a gene profile.
//' @return A numeric value.
//' @examples
//' ## alpha and beta subunits from the F-type ATP synthase.
//' data(fatp)
//' a <- t(fatp$atpPhylo['ATP5A1', ])
//' b <- t(fatp$atpPhylo['ATP5B', ])
//'
//' ## Person's correlation coefficient
//' corAB <- SimCor(a, b)
//' ## Jaccard similarity
//' jacAB <- SimJaccard(a, b)
//' ## Mutual information
//' MIAB <- SimMIBin(a, b)
//' MIABConti <- SimMIConti(a, b, bin = 10)
//' ## Hamming distance
//' hamAB <- DistHamming(a, b)
//' ## Eulidean distance
//' euAB <- DistEuclidean(a, b)
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @rdname simdist
//' @seealso SimDistBatch
//' @export
// [[Rcpp::export]]
double SimCor(arma::vec f,
              arma::vec t) {

  mat corMat = cor(f, t);
  double corVal = corMat(0, 0);

  return corVal;
}

//' @inheritParams SimCor
//' @rdname simdist
//' @export
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
//' @export
// [[Rcpp::export]]
arma::uword DistHamming(arma::vec f,
                        arma::vec t) {

  uword ham = sum(f != t);

  return ham;

}

//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double DistEuclidean(arma::vec f,
                     arma::vec t) {

  vec neq = f - t;

  double eu = sqrt(sum(square(neq)));

  return eu;
}

