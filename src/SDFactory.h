#ifndef SDFACTORY_H_
#define SDFACTORY_H_

#include <bigmemory/MatrixAccessor.hpp>
#include <memory>
#include "SDmeasure.h"

// [[Rcpp::depends(bigmemory)]]

//==============================
// SD Factory
//==============================
class SDFactory {
private:
  // p: The phylogenetic profile
  // idx: The index matrix
  // idxBig: The index big matrix
  const arma::mat p;
  const arma::umat idx;
  const arma::Mat<int> idxbig;
public:
  explicit SDFactory(const arma::mat& p,
                     const arma::Mat<int>& idxbig)
    : p(p), idxbig(idxbig) {};

  explicit SDFactory(const arma::mat& p,
                     const arma::umat& idx)
    : p(p), idx(idx) {};

  std::shared_ptr<SDmeasure> createSDFunc(Rcpp::List& attrs,
                                          Rcpp::List& arguments);
};


#endif
