#ifndef _DOLLO_H_
#define _DOLLO_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec MergeList(Rcpp::List x);

arma::uword CountRepeatIdx(arma::uword idx,
                           arma::uvec y);

arma::uvec InferGainNodes(Rcpp::List gainList);

arma::imat InferEdge(arma::umat edgeMat,
                     Rcpp::List tipPath,
                     arma::uvec pr);

arma::uword DolloDist(arma::umat edgeMat,
                      Rcpp::List tipPath,
                      arma::uvec pr1,
                      arma::uvec pr2);

#endif
