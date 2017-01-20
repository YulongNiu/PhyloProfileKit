#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "MI.h"

using namespace Rcpp;
using namespace arma;

//' @inheritParams SimCor
//' @param bin Integer.
//' @rdname simdist
//' @export
// [[Rcpp::export]]
double SimMIConti(arma::mat pairProfile,
                  arma::uword bin) {

  double MI;

  vec f = pairProfile.col(0);
  vec t = pairProfile.col(1);

  double n = f.n_elem;

  MI = Info(hist(f, bin), n) + Info(hist(t, bin), n) - Info(HistTwo(f, t, bin), n);

  return MI;

}


// [[Rcpp::export]]
double Info(arma::uvec x,
            double n) {

  uvec xPos = x(find(x != 0));

  vec xPosV = conv_to<vec>::from(xPos);

  vec p = xPosV / n;

  return -sum(p % log(p));
}

// [[Rcpp::export]]
arma::uvec HistTwo(arma::vec x,
                   arma::vec y,
                   arma::uword bin) {

  uvec count;

  uvec xIdx = FindInter(x, gInter(x, bin));
  uvec yIdx = FindInter(y, gInter(y, bin));

  // unique indices  elements
  uvec uYIdx = find_unique(yIdx, false);

  for (uword i = 0; i < uYIdx.n_elem; ++i) {
    uvec xV = xIdx(find(yIdx == yIdx(uYIdx(i))));
    count = join_cols(count, CountRepeat(xV));
  }

  return count;

}


// [[Rcpp::export]]
arma::uvec FindInter(arma::vec x,
                     arma::vec interval) {

  uword xN = x.n_elem;
  uvec idx(xN);

  for (uword i = 0; i < xN; ++i) {
    idx(i) = FindInterSingle(x(i), interval);
  }

  return idx;
}

// [[Rcpp::export]]
arma::uword FindInterSingle(double value,
                            arma::vec interval) {

  uword interN = interval.n_elem;

  for (uword i = 1; i < interN - 1; ++i) {
    if (value <= interval(i)) {
      return i;
    } else {}
  }

  return interN - 1;

}

// [[Rcpp::export]]
arma::vec gInter (arma::vec x,
                  arma::uword bin) {

  double minVal = x.min();
  double maxVal = x.max();
  double d = maxVal - minVal;
  vec b;

  if (d == 0) {
    d = abs(minVal);
    b = linspace<vec>(minVal - d / 1000, maxVal + d / 1000, bin + 1);
  } else {
    b = linspace<vec>(minVal, maxVal, bin + 1);
    b(0) -= d / 1000;
    b(bin) += d / 1000;
  }

  return b;

}


// [[Rcpp::export]]
arma::uvec CountRepeat (arma::uvec x) {

  vec u(1);
  uvec uNum(1);
  uword xL = x.n_elem;
  uword j;

  u(0) = x(0);
  uNum(0) = 1;

  for(uword i = 1; i < xL; ++i) {

    uword n = u.n_elem;

    for (j = 0; j < n; ++j) {
      if (x(i) == u(j)) {
        uNum(j) += 1;
        break;
      } else {}
    }

    if (j == n) {
      // add new unique elements
      u.resize(n + 1);
      u(n) = x(i);
      uNum.resize(n + 1);
      uNum(n) = 1;
    } else {}

  }

  return uNum;
}
