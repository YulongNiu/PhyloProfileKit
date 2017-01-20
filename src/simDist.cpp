#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simDist.h"

using namespace Rcpp;
using namespace arma;

//' Similarity or distance of paired phylogenetic profile
//'
//' SimCor(): Person's correlation coefficient.
//' SimJaccard(): Jaccard similarity.
//' SimMIBin(): Mutual information for binning data.
//' SimMIConti(): Mutual information for continuous data.
//' DistHamming(): Hamming distance.
//' DistEuclidean(): Euclidean distance.
//'
//' @title similarity and distance
//' @param pairProfile A paired phylogenetic profile, columns are genes and rows are species.
//' @return A numeric value.
//' @examples
//' ## alpha and beta subunits from the F-type ATP synthase.
//' data(fatp)
//' ab <- t(fatp$atpPhylo[c('ATP5A1', 'ATP5B'), ])
//'
//' ## Person's correlation coefficient
//' corAB <- SimCor(ab)
//' ## Jaccard similarity
//' jacAB <- SimJaccard(ab)
//' ## Mutual information
//' MIAB <- SimMIBin(ab)
//' MIABConti <- SimMIConti(ab, bin = 10)
//' ## Hamming distance
//' hamAB <- DistHamming(ab)
//' ## Eulidean distance
//' euAB <- DistEuclidean(ab)
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @rdname simdist
//' @seealso SimDistBatch
//' @export
// [[Rcpp::export]]
double SimCor(arma::mat pairProfile) {

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);

  mat corMat = cor(f, t);
  double corVal = corMat(0, 0);

  return corVal;
}

//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double SimJaccard(arma::mat pairProfile) {

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);
  vec combVec = f + 2*t;

  double A = sum(combVec == 3);

  double jac = A / (sum(f) + sum(t) - A);

  return jac;
}





//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
arma::uword DistHamming(arma::mat pairProfile) {

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);

  uword ham = sum(f != t);

  return ham;

}

//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double DistEuclidean(arma::mat pairProfile) {

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);

  vec neq = f - t;

  double eu = sqrt(sum(square(neq)));

  return eu;
}



