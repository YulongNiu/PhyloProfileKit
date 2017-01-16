#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "normSVD.h"

using namespace Rcpp;
using namespace arma;


//' SVD normalization
//'
//' Algorithm:
//'
//' Step1: rawBitM < bitCutoff to bitReset;
//'
//' Step2: In each row (species), x/max(x);
//'
//' Step3: L^2 SVD normalization.
//'
//' The core SVD normalization is retrieved from the SVD-Phy package with performance modification.
//' @title Singular value decomposition normalization of bit score matrix
//' @param rawBitM Raw bit score matrix.
//' @param bitCutoff Minimum value of the bit score.
//' @return
//'
//' SVDNorm(): SVD normalized bit score matrix.
//'
//' SVDPhy(): A L^2 normalized unitary matrix.
//'
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @inheritParams SVDPhy
//' @references \href{https://bitbucket.org/andrea/svd-phy}{SVD-Phy package}
//' @references \href{http://bioinformatics.oxfordjournals.org/content/suppl/2015/11/25/btv696.DC1/SVD-Phy-supplementary-material.docx}{SVD description}
//' @rdname normSVD
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix SVDNorm(Rcpp::NumericMatrix rawBitM,
                            double bitCutoff,
                            double bitReset,
                            double minConserve,
                            double trimming) {

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
                                 bitReset = bitReset,
                                 minConserve = minConserve,
                                 trimming = trimming);

  return resultM;
}



//' @param bitM Bit score matrix, for example the BLASTP or STRING bit scores. It is a named numeric matrix, columns are species and rows are genes.
//' @param bitReset Reset the bit score for ones lower than the `bitCutoff`.
//' @param trimming A percentages top unitary matrix.
//' @param minConserve Minimum number of homologous. The proteins with homologous less than this value are discarded.
//' @rdname normSVD
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix SVDPhy(Rcpp::NumericMatrix bitM,
                           double bitReset,
                           double minConserve,
                           double trimming) {

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

  // filter genes (rows)
  uvec hasV = sum(bitMArma > bitReset, 1);
  uvec filteredIdx = find(hasV >= minConserve);
  U = U.rows(filteredIdx);

  // trim species
  mat UTrimmed = U.cols(regspace<uvec>(0, round(trimming * U.n_cols) - 1));

  // L^2 normalization
  mat normU = normalise(UTrimmed, 2, 1);

  // rownames
  NumericVector filteredIdxCpp(filteredIdx.begin(), filteredIdx.end());
  NumericMatrix resultM(normU.n_rows, normU.n_cols, normU.begin());
  rownames(resultM) = rn[filteredIdxCpp];

  return resultM;
}

