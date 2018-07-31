#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <bigmemory/MatrixAccessor.hpp>

#include <algorithm>

#include "SDFactory.h"
#include "SDmeasure.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo, RcppParallel, bigmemory)]]

template <typename T>
struct BatchCore : public Worker {

  const mat& p;
  const T& idx;
  std::shared_ptr<SDmeasure>& sdfunc;
  NumericVector& res;

  BatchCore(const arma::mat& p,
            const T& idx,
            std::shared_ptr<SDmeasure>& sdfunc,
            Rcpp::NumericVector& res)
    : p(p), idx(idx), sdfunc(sdfunc), res(res) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; ++i) {

      uword fidx = idx(i, 0);
      uword tidx = idx(i, 1);

      rowvec frow = p.row(fidx-1);
      vec f(frow.begin(), frow.n_elem);
      rowvec trow = p.row(tidx-1);
      vec t(trow.begin(), trow.n_elem);

      res[i] = sdfunc->calcSD(f, t);
    }
  }
};



// [[Rcpp::export]]
Rcpp::NumericVector BatchMat(const arma::mat p,
                             const arma::umat idx,
                             Rcpp::List attrs,
                             Rcpp::List arguments) {

  unsigned long n = idx.n_rows;

  // allocate the result vector we will return
  NumericVector res(n);

  // create the worker
  std::shared_ptr<SDmeasure> sdfunc = SDFactory(p, idx).createSDFunc(attrs, arguments);
  BatchCore<arma::umat>* batchWorker = new BatchCore<arma::umat>(p, idx, sdfunc, res);

  // call it with parallelFor
  parallelFor(0, n, (*batchWorker));
  delete batchWorker;
  batchWorker = NULL;

  return res;
}



// [[Rcpp::export]]
Rcpp::NumericVector BatchBigmat(const arma::mat p,
                                SEXP idx,
                                Rcpp::List attrs,
                                Rcpp::List arguments) {

  XPtr<BigMatrix> bigidx_(idx);
  unsigned long n = bigidx_->nrow();
  Mat<int> bigidx = Mat<int>((int *)bigidx_->matrix(), bigidx_->nrow(), bigidx_->ncol(), false);

  // allocate the result vector we will return
  NumericVector res(n);

  // create the worker
  std::shared_ptr<SDmeasure> sdfunc = SDFactory(p, bigidx).createSDFunc(attrs, arguments);
  BatchCore<arma::Mat<int>>* batchWorker = new BatchCore<arma::Mat<int>>(p, bigidx, sdfunc, res);

  // call it with parallelFor
  parallelFor(0, n, (*batchWorker));
  delete batchWorker;
  batchWorker = NULL;

  return res;
}
