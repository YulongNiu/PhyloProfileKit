#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "phylo2Mat.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec Phylo2Mat(arma::umat edgeMat,
                    arma::uword tipNum) {

  // construct init x location
  umat initE = InitEdge(edgeMat, tipNum);
  vec initX = InitX(initE, tipNum);

  return x;
}


// [[Rcpp::export]]
arma::vec InitX(arma::umat m,
                arma::uword tipNum) {
  uword nr = m.n_rows;
  vec xloc = zeros<vec>(nr);

  for (uword i = 0; i < nr; ++i) {
    if (m(i, 1) < tipNum) {
      xloc(i) = m(i, 1);
    } else {}
  }

  return xloc;
}

// [[Rcpp::export]]
arma::umat InitEdge(arma::umat edgeMat,
                    arma::uword tipNum) {
  urowvec root(2);
  root.fill(tipNum + 1);

  edgeMat.insert_rows(edgeMat.n_rows, root);

  return edgeMat;
}
