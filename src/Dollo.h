#ifndef _DOLLO_H_
#define _DOLLO_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec MergeList(Rcpp::List x);

arma::uvec InferGainNodes(Rcpp::List gainList);

arma::imat InferEdge(arma::umat edgeMat,
                     Rcpp::List tipPath,
                     Rcpp::NumericVector pr);

arma::uword DolloDist(arma::umat edgeMat,
                      Rcpp::List tipPath,
                      Rcpp::NumericVector pr1,
                      Rcpp::NumericVector pr2);

#endif
