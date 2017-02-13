#ifndef _PHYLO2MAT_H_
#define _PHYLO2MAT_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::vec Phylo2Mat(arma::umat edgeMat,
                    arma::uword tipNum);

arma::vec InitX(arma::umat m,
                arma::uword tipNum);

arma::umat InitEdge(arma::umat edgeMat,
                    arma::uword tipNum);

#endif
