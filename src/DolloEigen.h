#ifndef _DOLLOEIGEN_H_
#define _DOLLOEIGEN_H_

#include <RcppEigen.h>
#include <Rcpp.h>

int CountRepeatIdxEigen(int idx,
                        Eigen::VectorXi y);

void TestLen(Eigen::VectorXi y);

#endif
