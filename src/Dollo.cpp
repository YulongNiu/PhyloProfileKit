#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "Dollo.h"

using namespace Rcpp;
using namespace arma;

//' Use Dollo's parsimony distance to evaluate similarity patterns.
//'
//' Algorithm:
//' 1. Find the latest common ancestors (LCA) for all gain tips.
//' 2. Select the nodes from each tip to the LCA.
//'
//' InferGainNodes(): Infer the present nodes and tips.
//'
//' DolloDist(): Infer the Dollo's parsimony distance.
//'
//' @title Dollo's parsimony distance
//' @param gainList A list of ancestors of each gain tips. In each elements, the first one is the root.
//' @return
//'
//' InferGainNodes(): A vector (with tips) of 1 and 0.
//'
//' DolloDist(): An integer.
//'
//' @examples
//' ## example tree
//' library('ape')
//' testTreeText <- '((((t1, t2),(t3, t4)),t5), (t6, (t7, t8)));'
//' testTree <- read.tree(text = testTreeText)
//' plot(testTree)
//'
//' ## node path
//' pathList <- nodepath(testTree)
//' InferGainNodes(pathList[1])
//' InferGainNodes(pathList[1:3])
//'
//' InferEdge(testTree$edge, pathList, c(0, 0, 0, 0, 1, 1, 1, 1))
//'
//' DolloDist(testTree$edge, pathList, c(0, 0, 0, 0, 1, 1, 1, 1), c(0, 0, 0, 0, 0, 1, 1, 1))
//' DolloDist(testTree$edge, pathList, c(1, 1, 1, 1, 1, 0, 0, 1), c(0, 0, 0, 0, 1, 0, 0, 1))
//' @author Yulong Niu \email{niuylscu@@gmail.com}
//' @rdname dollo
//' @export
// [[Rcpp::export]]
arma::uvec InferGainNodes(Rcpp::List gainList) {

  uvec gainVec;

  uword gainNum = gainList.size();
  uvec gainA = MergeList(gainList);
  uword gainANum = gainA.size();

  uword root = gainA[0];
  uword ancestor = 0;

  if (gainNum > 1) {
    for (uword i = 0; i < gainANum; ++i) {
      uword eachRepeatNum = CountRepeatIdx(i, gainA);
      if (eachRepeatNum < gainNum) {
        ancestor = gainA[i - 1];
        break;
      } else {}
    }

    gainA = unique(gainA);
    gainVec = gainA.elem(find(gainA < root || gainA >= ancestor));

  } else {
    gainVec = gainA.elem(find(gainA < root));
  }

  return gainVec;
}


//' @export
// [[Rcpp::export]]
arma::uvec InferGainNodes2(Rcpp::List gainList) {

  uvec gainVec;

  uword gainNum = gainList.size();
  uvec path1st = gainList[0];
  uvec gainA = MergeList(gainList);
  uword gainANum = gainA.size();

  uword root = gainA[0];
  uword ancestor = 0;

  if (gainNum > 1) {

    uvec gainAUni = unique(gainA);
    uword nodeNum = gainAUni.size();
    uvec (nodeNum, fill::zeros);

    for (uowrd i = 0; i < nodeNum; ++i) {
      
    }
  } else {
    gainVec = gainA.elem(find(gainA < root));
  }

  if (gainNum > 1) {
    for (uword i = 0; i < gainANum; ++i) {
      uword eachRepeatNum = CountRepeatIdx(i, gainA);
      if (eachRepeatNum < gainNum) {
        ancestor = gainA[i - 1];
        break;
      } else {}
    }

    gainA = unique(gainA);
    gainVec = gainA.elem(find(gainA < root || gainA >= ancestor));

  } else {
    gainVec = gainA.elem(find(gainA < root));
  }

  return gainVec;
}


//' @param edgeMat A edge mat could be generated from the "ape" package.The first row should be (root --> nodes).
//' @param tipPath  A list of ancestors of each gain tips. In each elements, the first one is the root.
//' @param pr A numeric vector indicates the "presence-absence" pattern."pr" should be in the same order with the tips of tree
//' @return A numeric edge present and absent matrix.
//' @rdname dollo
//' @export
// [[Rcpp::export]]
arma::imat InferEdge(arma::umat edgeMat,
                     Rcpp::List tipPath,
                     arma::uvec pr) {

  // gain and loss matrix
  imat glMat(edgeMat.n_rows, 2, fill::zeros);
  uword gainTipNum = accu(pr);

  if (gainTipNum > 0 && gainTipNum < pr.size()) {
    // gain tips
    List gainTips = tipPath[as<IntegerVector>(wrap(find(pr == 1)))];

    // infer gain nodes and tips
    uvec gains = InferGainNodes(gainTips);

    for (uword i = 0; i < gains.size(); ++i) {
      glMat.elem(find(edgeMat == gains[i])).ones();
    }
  }
  else if (gainTipNum == pr.size()) {
    // all gain tips
    glMat.ones();
  }
  // if no gain tips, then all zero
  else{}

  return glMat;
}


//' @param pr1 1st "presence-absence" pattern.
//' @param pr2 2nd "presence-absence" pattern.
//' @inheritParams InferEdge
//' @return A number indicating the Dollo's parsimony distance.
//' @rdname dollo
//' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/?term=17535793}
//' @export
// [[Rcpp::export]]
arma::uword DolloDist(arma::umat edgeMat,
                      Rcpp::List tipPath,
                      arma::uvec pr1,
                      arma::uvec pr2) {

  imat glMat1 = InferEdge(edgeMat, tipPath, pr1);
  imat glMat2 = InferEdge(edgeMat, tipPath, pr2);

  uword dist = accu(abs((glMat1.col(0) - glMat1.col(1)) - (glMat2.col(0) - glMat2.col(1))));

  return dist;
}


//' @param x A list only contains unsigned interger vectors.
//' @return A merged numeric vector.
//' @rdname dollo
//' @keywords internal
// [[Rcpp::export]]
arma::uvec MergeList(Rcpp::List x) {

  uword len = x.size();
  uvec mergedX;

  for (uword i = 0; i < len; ++i) {
    uvec eachX = x[i];
    mergedX = join_cols(mergedX, eachX);
  }

  return mergedX;

}


//' @param idx Index
//' @param y A unsigned integer vector
//' @return Repeat number of y[idx] in y.
//' @rdname dollo
//' @keywords internal
// [[Rcpp::export]]
arma::uword CountRepeatIdx(arma::uword idx,
                           arma::uvec y) {

  uword repeatNum = 0;
  uword yNum = y.size();
  uword target = y[idx];

  for (uword i = 0; i < yNum; ++i) {
    if (target == y[i]) {
      ++repeatNum;
    } else {}
  }

  return repeatNum;
}
