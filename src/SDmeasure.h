#ifndef SDMEASURE_H_
#define SDMEASURE_H_

#include <RcppArmadillo.h>

typedef double (*funcPtr)(const arma::vec &f, const arma::vec &t);

using namespace std;
using namespace arma;

class SDmeasure {
public:
  virtual ~SDmeasure() {};
  virtual double calcSD(const arma::vec &f, const arma::vec &t) = 0;
};

#endif

