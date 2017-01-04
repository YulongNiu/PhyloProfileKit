#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simDistCpp.h"

using namespace Rcpp;
using namespace arma;


//' Similarity or distance of paired phylogenetic profile
//'
//' SimCor(): Person's correlation coefficient.
//' SimJaccard(): Jaccard similarity.
//' SimMI(): Mutual information.
//' SimMIConti(): Mutual information for continuous variables
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
//' MIAB <- SimMI(ab)
//' MIABConti <- SimMIConti(ab)
//' ## Hamming distance
//' hamAB <- DistHamming(ab)
//' ## Eulidean distance
//' euAB <- DistEuclidean(ab)
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @importFrom stats cor
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

  uvec inter = find(combVec == 3);

  double A = inter.n_elem;

  double jac = A / (sum(f) + sum(t) - A);

  return jac;
}


//' @inheritParams SimCor
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double SimMI(arma::mat pairProfile) {

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);
  vec combVec = f + 2*t;

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
arma::uword DistHamming(arma::mat pairProfile) {

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);

  uvec neq = find(f != t);

  uword ham = neq.n_elem;

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


// // [[Rcpp::export]]
// arma::umat testHist(arma::mat tmp1) {
//   return hist(tmp1, linspace<vec>(tmp1.min(), tmp1.max(), 10));
// }
