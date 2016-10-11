#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "DolloEigen.h"

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
// // [[Rcpp::export]]
// int CountRepeatIdxEigen(int idx,
//                         Eigen::VectorXi y){

//   int repeatNum = 0;
//   int yNum = y.size();
//   int target = y(idx);

//   for (int i = 0; i < yNum; ++i) {
//     if (target == y(i)) {
//       ++repeatNum;
//     } else {}
//   }

//   return repeatNum;
// }



// [[Rcpp::export]]
void TestLen(VectorXd y){

  for(int i = 0; i < 10; ++i) {
    std::cout << y(i) << std::endl;
  }

}

// [[Rcpp::export]]
VectorXd getEigenValues(Map<MatrixXd> M) {
  SelfAdjointEigenSolver<MatrixXd> es(M);
  return es.eigenvalues();
}
