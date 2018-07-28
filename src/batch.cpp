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

struct BatchCoreMat : public Worker {

  const mat& p;
  const umat& idx;
  std::shared_ptr<SDmeasure>& sdfunc;
  vec& res;

  BatchCoreMat(const arma::mat& p,
               const arma::umat& idx,
               std::shared_ptr<SDmeasure>& sdfunc,
               arma::vec& res)
    : p(p), idx(idx), sdfunc(sdfunc), res(res) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; ++i) {

      uword fidx = idx(i, 0);
      uword tidx = idx(i, 1);

      rowvec f = p.row(fidx-1);
      rowvec t = p.row(tidx-1);

      res(i) = sdfunc->calcSD(f, t);
    }
  }
};


struct BatchCoreBigmat : public Worker {

  const mat& p;
  SEXP idx;
  std::shared_ptr<SDmeasure>& sdfunc;
  vec& res;

  BatchCoreBigmat(const arma::mat& p,
                  SEXP idx,
                  std::shared_ptr<SDmeasure>& sdfunc,
                  arma::vec& res)
    : p(p), idx(idx), sdfunc(sdfunc), res(res) {}

  void operator()(std::size_t begin, std::size_t end) {

    XPtr<BigMatrix> bigidx(idx);
    Mat<int> bigidxarma = Mat<int>((int *)bigidx->matrix(), bigidx->nrow(), bigidx->ncol(), false);

    for (std::size_t i = begin; i < end; ++i) {

      uword fidx = bigidxarma(i, 0);
      uword tidx = bigidxarma(i, 1);

      rowvec f = p.row(fidx-1);
      rowvec t = p.row(tidx-1);

      res(i) = sdfunc->calcSD(f, t);
    }
  }
};



// [[Rcpp::export]]
arma::vec BatchMat(const arma::mat p,
                   const arma::umat idx,
                   Rcpp::List attrs,
                   Rcpp::List arguments) {

  unsigned long n = idx.n_rows;

  // allocate the result vector we will return
  vec res(n);

  // create the worker
  std::shared_ptr<SDmeasure> sdfunc = SDFactory(p, idx).createSDFunc(attrs, arguments);
  BatchCoreMat* batchWorker = new BatchCoreMat(p, idx, sdfunc, res);

  // call it with parallelFor
  parallelFor(0, n, (*batchWorker));
  delete batchWorker;
  batchWorker = NULL;

  return res;
}



// [[Rcpp::export]]
arma::vec BatchBigmat(const arma::mat p,
                      SEXP idx,
                      Rcpp::List attrs,
                      Rcpp::List arguments) {

  XPtr<BigMatrix> bigidx(idx);
  unsigned long n = bigidx->nrow();

  // allocate the result vector we will return
  vec res(n);

  // create the worker
  std::shared_ptr<SDmeasure> sdfunc = SDFactory(p, idx).createSDFunc(attrs, arguments);
  BatchCoreBigmat* batchWorker = new BatchCoreBigmat(p, idx, sdfunc, res);

  // call it with parallelFor
  parallelFor(0, n, (*batchWorker));
  delete batchWorker;
  batchWorker = NULL;

  return res;
}
