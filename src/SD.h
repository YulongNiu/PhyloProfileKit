#ifndef SD_H_
#define SD_H_

#include "SDmeasure.h"

class SimJaccard : public SDmeasure {
public:
  double calcSD(const arma::rowvec& f,
                const arma::rowvec& t) {

    vec combVec = f + 2*t;
    double A = sum(combVec == 3);

    return  A / (sum(f) + sum(t) - A);
  }
};


class DistHamming : public SDmeasure {
public:
  double calcSD(const arma::rowvec& f,
                const arma::rowvec& t) {

    return sum(f != t);
  }
};


class DistManhattan : public SDmeasure {
public:
  double calcSD(const arma::rowvec& f,
                const arma::rowvec& t) {

    return sum(abs(f - t));
  }
};


class DistEuclidean : public SDmeasure {
public:
  double calcSD(const arma::rowvec& f,
                const arma::rowvec& t) {

    return sqrt(sum(square(f - t)));
  }
};


class DistMinkowski : public SDmeasure {
private:
  unsigned int p;
public:
  explicit DistMinkowski (unsigned int p) {
    this->p = p;
  }
  ~DistMinkowski() {}
  double calcSD(const arma::rowvec& f,
                const arma::rowvec& t) {
    return pow(arma::accu(arma::pow(arma::abs(f - t), this->p)), 1.0 / this->p);
  }
};


#endif
