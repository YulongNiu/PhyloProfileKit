#ifndef SD_H_
#define SD_H_

#include "SDmeasure.h"
#include "MI.h"


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


//' @param p \code{numeric} value indicating p-norm parameter.
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
//======================================
//  Binary mutual information (similarity)
//======================================
class SimMIBin : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {

    vec combVec = f + 2*t;

    double N = f.n_elem;
    double A = sum(combVec == 3);
    double B = sum(combVec == 1);
    double C = sum(combVec == 2);
    double D = N - A - B - C;

    double I = eachMI(A, B, C, N) + eachMI(B, A, D, N) + eachMI(C, A, D, N) + eachMI(D, C, B, N);

    return I;
  }
};

//' @param bin A positive \code{Integer} indicating the bin.
//' @inheritParams SimCor
//' @rdname simdist
//' @keywords internal
//=============================================
//  Continuous mutual information (similarity)
//=============================================
class SimMIConti : public SDmeasure {
private:
  uword bin;
public:
  explicit SimMIConti (arma::uword bin) {
    this->bin = bin;
  }
  ~SimMIConti() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {

    double n = f.n_elem;
    double MI = Info(hist(f, this->bin), n) + Info(hist(t, this->bin), n) - Info(HistTwo(f, t, this->bin), n);

    return MI;
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


// class SDCustom : public SDmeasure {
// private:
//   funcPtr func;
//   Rcpp::List funcargs;
// public:
//   explicit SDCustom (funcPtr function, Rcpp::List funcargs) :
//     func(function), funcargs(funcargs) {
//     this->func = function;
//     this->funcargs = funcargs;
//   };
//   ~SDCustom () {}
//   double calcSD(const arma::vec &f,
//                 const arma::vec &t) {
//     return func(f, t, funcargs);
//   }
// };

#endif
