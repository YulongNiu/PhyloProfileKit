#ifndef SDMEASURE_H_
#define SDMEASURE_H_

#include <RcppArmadillo.h>

typedef double (*funcPtr)(const arma::rowvec &f, const arma::rowvec &t);

using namespace std;
using namespace arma;

class SDmeasure {
public:
  virtual ~SDmeasure() {};
  virtual double calcSD(const arma::rowvec &f, const arma::rowvec &t) = 0;
};

#endif

