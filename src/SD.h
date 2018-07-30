#ifndef SD_H_
#define SD_H_

#include "SDmeasure.h"


//==============================================
// Person's correlation coefficient (similarity)
//==============================================
class SimCor : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat corMat = cor(f, t);
    return corMat(0, 0);
  }
};


//=======================
// Jaccard similarity
//=======================
class SimJaccard : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    vec combVec = f + 2*t;
    double A = sum(combVec == 3);
    return  A / (sum(f) + sum(t) - A);
  }
};

//=======================
// Hamming distance
//=======================
class DistHamming : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return sum(f != t);
  }
};


//=======================
// Manhattan distance
//=======================
class DistManhattan : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return sum(abs(f - t));
  }
};


//=======================
// Euclidean distance
//=======================
class DistEuclidean : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return sqrt(sum(square(f - t)));
  }
};


//=======================
// Minkowski distance
//=======================
class DistMinkowski : public SDmeasure {
private:
  unsigned int p;
public:
  explicit DistMinkowski (unsigned int p) {
    this->p = p;
  }
  ~DistMinkowski() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return pow(arma::accu(arma::pow(arma::abs(f - t), this->p)), 1.0 / this->p);
  }
};


//============================
// Custom similarity/distance
//============================
class SDCustom : public SDmeasure {
private:
  funcPtr func;
public:
  explicit SDCustom (funcPtr function) : func(function) {
    this->func = function;
  };
  ~SDCustom () {}
  double calcSD(const arma::vec &f,
                const arma::vec &t) {
    return func(f, t);
  }
};

#endif
