#ifndef _COLLAPSETREE_H_
#define _COLLAPSETREE_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat CollapseTree(arma::mat edgeMat,
                       arma::uword tipNum,
                       arma::vec f,
                       arma::vec t);

bool isTwoRowsEqual(arma::mat m);

#endif
