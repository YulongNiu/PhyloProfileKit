#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "normSVD.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix SVDPhy(Rcpp::NumericMatrix bitM,
                           double trimming,
                           double minConserve) {

  //keep rownames
  CharacterVector rn = rownames(bitM);
  mat bitMArma(bitM.begin(), bitM.nrow(), bitM.ncol(), false);

  uvec dimV(2);
  dimV(0) = bitMArma.n_rows;
  dimV(1) = bitMArma.n_cols;

  // SVD
  mat U;
  vec s;
  mat V;
  svd(U, s, V, bitMArma);
  U = U.cols(0, dimV.min() - 1);

  // filter genes
  uvec hasV(bitMArma.n_rows);
  for (uword i = 0; i < bitMArma.n_rows; ++i) {
    uvec counterIdx  = find(bitMArma.row(i) > 0);
    hasV(i) = counterIdx.n_elem;
  }
  uvec filteredIdx = find(hasV > minConserve);
  U = U.rows(filteredIdx);

  // trim species
  mat UTrimmed = U.cols(regspace<uvec>(0, round(trimming * U.n_cols) - 1));

  // L^2 normalization
  UTrimmed.each_row([](rowvec& r){
      r /= sqrt(sum(square(r)));
    });

  // rownames
  NumericVector filteredIdxCpp(filteredIdx.begin(), filteredIdx.end());
  NumericMatrix resultM(UTrimmed.n_rows, UTrimmed.n_cols, UTrimmed.begin());
  rownames(resultM) = rn[filteredIdxCpp];

  return resultM;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix SVDNorm(Rcpp::NumericMatrix rawBitM,
                           double bitCutoff,
                           double bitReset,
                           double trimming,
                           double minConserve) {

  // keep rownames
  CharacterVector rn = rownames(rawBitM);
  mat rawBitMArma(rawBitM.begin(), rawBitM.nrow(), rawBitM.ncol(), false);

  // step1: rawBitM < bitCutoff to bitReset
  rawBitMArma.for_each([bitCutoff, bitReset](mat::elem_type& val) {
      val = val < bitCutoff ? bitReset : val;
    });

  // step2: In each row (species), x/max(x)
  rawBitMArma.each_row([](rowvec& r) {
      r /= r.max();
    });

  // step3: L^2 SVD normalization
  NumericMatrix norP(rawBitMArma.n_rows, rawBitMArma.n_cols, rawBitMArma.begin());
  rownames(norP) = rn;

  NumericMatrix resultM = SVDPhy(norP,
                                 trimming = trimming,
                                 minConserve = minConserve);

  return resultM;
}

