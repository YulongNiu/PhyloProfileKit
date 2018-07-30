#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "MI.h"

using namespace Rcpp;
using namespace arma;


//' Utilities for MI
//'
//' \code{eachMI()}: Info for a cell.
//' \code{Info()}: Entropy.
//' \code{HistTwo()}: Joint counts of two vectors.
//' \code{FindInter()}: Interval indices of a vector.
//' \code{FindInterSingle()}: Interval index of a value.
//' \code{gInter()}: Generate an interval vector, inspired from the \code{cut()} function of the \code{base} package.
//' \code{CountRepeat()}: Repeat counts of a vector.
//'
//' @param p1, p2, p3: Counts of variables in cells.
//' @param n Total variables.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname utilities-MI
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


//' @param v Histogram of counts.
//' @inheritParams eachMI
//' @rdname utilities-MI
//' @keywords internal
// [[Rcpp::export]]
double Info(arma::uvec v,
            double n) {

  uvec xPos = v(find(v != 0));

  vec xPosV = conv_to<vec>::from(xPos);

  vec p = xPosV / n;

  return -sum(p % log(p));
}

//' @param bin A positive \code{integer} indicating the bin.
//' @param x, y \code{numeric vector}.
//' @rdname utilities-MI
//' @keywords internal
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

//' @inheritParams HistTwo
//' @param internal Interval numeric vector.
//' @rdname utilities-MI
//' @keywords internal
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

//' @inheritParams FindInter
//' @param value Number.
//' @rdname utilities-MI
//' @keywords internal
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

//' @inheritParams HistTwo
//' @rdname utilities-MI
//' @keywords internal
// [[Rcpp::export]]
arma::vec gInter(arma::vec x,
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

//' @inheritParams HistTwo
//' @rdname utilities-MI
//' @keywords internal
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
