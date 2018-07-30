#ifndef SD_H_
#define SD_H_

#include "SDmeasure.h"


//' Similarity or distance of paired phylogenetic profile
//'
//' \code{SimCor()}: Person's correlation coefficient.
//'
//' \code{SimJaccard()}: Jaccard similarity.
//'
//' \code{SimMIBin()}: Mutual information for binning data.
//'
//' \code{SimMIConti()}: Mutual information for continuous data.
//'
//' \code{DistHamming()}: Hamming distance.
//'
//' \code{DistEuclidean()}: Euclidean distance.
//'
//' @title Similarity and distance
//' @param f Numeric vector indicating a gene profile.
//' @param t Numeric vector indicating a gene profile.
//' @return A numeric value.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname simdist
//' @keywords internal
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


//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
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

//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
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


//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
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


//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
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


//' @param p The p-norm parameter.
//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
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

//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
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
