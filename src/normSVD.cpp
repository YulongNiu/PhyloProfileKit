#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "normSVD.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat SVDArma(arma::mat bitM,
                  double trimming,
                  double minConserve) {

  uvec dimV(2);
  dimV(0) = bitM.n_rows;
  dimV(1) = bitM.n_cols;

  // SVD
  mat U;
  vec s;
  mat V;
  svd(U, s, V, bitM);
  U = U.cols(0, dimV.min() - 1);

  // filter genes
  vec filteredIdx(bitM.n_rows);
  for (uword i = 0; i < bitM.n_rows; ++i) {
    uvec counterIdx  = find(bitM.row(i) > 0);
    filteredIdx(i) = counterIdx.n_elem;
  }
  U = U.rows(find(filteredIdx > minConserve));

  // trim species
  uword trimNum = round(trimming * U.n_cols);
  mat UTrimmed = U.cols(linspace<uvec>(0, trimNum - 1, trimNum));

  // L^2 normalization
  mat resultM(size(UTrimmed));
  for (uword i = 0; i < UTrimmed.n_rows; ++i) {
    rowvec eachRow = UTrimmed.row(i);
    resultM.row(i) = eachRow/sqrt(sum(square(eachRow)));
  }

  return resultM;
}
