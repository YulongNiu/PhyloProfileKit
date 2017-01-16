#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "normNPP.h"

using namespace Rcpp;
using namespace arma;


//' NPP normalization
//'
//' Algorithm:
//'
//' Step1: rawBitM < hitCutoff to hitReset;
//'
//' Step2: filter genes without enough homologys
//'
//' Step3: in each row (species), log2(x/max(x));
//'
//' Step4: z-score for each column.
//' 
//' @title z-score normalization of phylogenetic profile
//' @return
//'
//' NPPNorm(): NPP normalized bit score matrix.
//'
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @inheritParams SVDNorm
//' @references \href{http://www.nature.com/nature/journal/v493/n7434/extref/nature11779-s1.pdf}{NPP description}
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix NPPNorm(Rcpp::NumericMatrix rawBitM,
                            double bitCutoff,
                            double bitReset,
                            double minConserve) {
  // keep names
  CharacterVector rn = rownames(rawBitM);
  CharacterVector cn = colnames(rawBitM);
  mat rawBitMArma(rawBitM.begin(), rawBitM.nrow(), rawBitM.ncol(), false);

  // step1: rawBitM < bitCutoff to bitReset
  rawBitMArma.for_each([bitCutoff, bitReset](mat::elem_type& val) {
      val = val < bitCutoff ? bitReset : val;
    });

  // step2: filter genes without enough homologys
  uvec hasV = sum(rawBitMArma > bitReset, 1);
  uvec filteredIdx = find(hasV >= minConserve);
  rawBitMArma = rawBitMArma.rows(filteredIdx);

  // step3: In each row (species), log2(x/max(x))
  rawBitMArma.each_row([](rowvec& r) {
      r = log2(r/r.max());
    });

  // step4: z-score for each column
  rawBitMArma.each_col([](vec& v) {
      v = (v - mean(v))/stddev(v);
    });

  // reassign names
  NumericVector filteredIdxCpp(filteredIdx.begin(), filteredIdx.end());
  NumericMatrix resultM(rawBitMArma.n_rows, rawBitMArma.n_cols, rawBitMArma.begin());
  rownames(resultM) = rn[filteredIdxCpp];
  colnames(resultM) = cn;

  return resultM;
}

