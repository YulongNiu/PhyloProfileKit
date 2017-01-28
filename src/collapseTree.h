#ifndef _COLLAPSETREE_H_
#define _COLLAPSETREE_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::umat CollapseTree(arma::umat edgeMat,
                        arma::uword tipNum,
                        arma::uvec f,
                        arma::uvec t);

bool isTwoRowsEqual(arma::umat m);

#endif
