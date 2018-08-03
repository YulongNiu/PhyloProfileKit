#ifndef SD_H_
#define SD_H_

#include "SDmeasure.h"
#include "MI.h"
#include "collapseTree.h"
#include "simDist.h"

//==============================================
// Person's correlation coefficient (similarity)
//==============================================

class SimCor : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return SimCor_(f, t);
  }
};

//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//=========================================================
// Collapsed Person's correlation coefficient (similarity)
//=========================================================
class SimCorCollapse : public SDmeasure {
private:
  mat edgeMat;
  uword tipNum;
public:
  explicit SimCorCollapse(arma::mat edgeMat, arma::uword tipNum) {
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~SimCorCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return SimCor_(fnew, tnew);
  }
};


class SimJaccard : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return  SimJaccard_(f, t);
  }
};


//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//==============================
// Collapsed Jaccard similarity
//==============================
class SimJaccardCollapse : public SDmeasure {
private:
  mat edgeMat;
  uword tipNum;
public:
  explicit SimJaccardCollapse(arma::mat edgeMat, arma::uword tipNum) {
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~SimJaccardCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return SimJaccard_(fnew, tnew);
  }
};


//==============================
// Hamming distance
//==============================
class DistHamming : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return DistHamming_(f, t);
  }
};


//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//==============================
// Collapsed Hamming distance
//==============================
class DistHammingCollapse : public SDmeasure {
private:
  mat edgeMat;
  uword tipNum;
public:
  explicit DistHammingCollapse(arma::mat edgeMat, arma::uword tipNum) {
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~DistHammingCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return DistHamming_(fnew, tnew);
  }
};



//=======================
// Manhattan distance
//=======================
class DistManhattan : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return DistManhattan_(f, t);
  }
};


//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//==============================
// Collapsed Manhattan distance
//==============================
class DistManhattanCollapse : public SDmeasure {
private:
  mat edgeMat;
  uword tipNum;
public:
  explicit DistManhattanCollapse(arma::mat edgeMat, arma::uword tipNum) {
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~DistManhattanCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return DistManhattan_(fnew, tnew);
  }
};

//=======================
// Euclidean distance
//=======================
class DistEuclidean : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return DistEuclidean_(f, t);
  }
};


//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//==============================
// Collapsed Euclidean distance
//==============================
class DistEuclideanCollapse : public SDmeasure {
private:
  mat edgeMat;
  uword tipNum;
public:
  explicit DistEuclideanCollapse(arma::mat edgeMat, arma::uword tipNum) {
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~DistEuclideanCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return DistEuclidean_(fnew, tnew);
  }
};


//=======================
// Minkowski distance
//=======================
class DistMinkowski : public SDmeasure {
private:
  uword p;
public:
  explicit DistMinkowski (uword p) {
    this->p = p;
  }
  ~DistMinkowski() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return DistMinkowski_(f, t, this->p);
  }
};

//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//==============================
// Collapsed Minkowski distance
//==============================
class DistMinkowskiCollapse : public SDmeasure {
private:
  uword p;
  mat edgeMat;
  uword tipNum;
public:
  explicit DistMinkowskiCollapse(uword p, arma::mat edgeMat, arma::uword tipNum) {
    this->p = p;
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~DistMinkowskiCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return DistMinkowski_(fnew, tnew, this->p);
  }
};



//=========================================
//  Binary mutual information (similarity)
//=========================================
class SimMIBin : public SDmeasure {
public:
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    return SimMIBin_(f, t);
  }
};


//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//==================================================
// Collapsed Binary mutual information (similarity)
//==================================================
class SimMIBinCollapse : public SDmeasure {
private:
  mat edgeMat;
  uword tipNum;
public:
  explicit SimMIBinCollapse(arma::mat edgeMat, arma::uword tipNum) {
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~SimMIBinCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return SimMIBin_(fnew, tnew);
  }
};


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
    return SimMIConti_(f, t, this->bin);
  }
};

//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//=====================================================
// Collapsed Continuous mutual information (similarity)
//=====================================================
class SimMIContiCollapse : public SDmeasure {
private:
  uword bin;
  mat edgeMat;
  uword tipNum;
public:
  explicit SimMIContiCollapse(uword bin, arma::mat edgeMat, arma::uword tipNum) {
    this->bin = bin;
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~SimMIContiCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return SimMIConti_(fnew, tnew, this->bin);
  }
};

//============================
// Custom similarity/distance
//============================
class SDCustom : public SDmeasure {
private:
  funcPtr func;
public:
  explicit SDCustom (funcPtr func) {
    this->func = func;
  }
  ~SDCustom () {}
  double calcSD(const arma::vec &f,
                const arma::vec &t) {
    return func(f, t);
  }
};


//' @inheritParams SimCor_
//' @rdname simdist
//' @keywords internal
//======================================
// Collapsed custom similarity/distance
//======================================
class SDCustomCollapse : public SDmeasure {
private:
  funcPtr func;
  mat edgeMat;
  uword tipNum;
public:
  explicit SDCustomCollapse(funcPtr func,
                            arma::mat edgeMat,
                            arma::uword tipNum) {
    this->func = func;
    this->edgeMat = edgeMat;
    this->tipNum = tipNum;
  }
  ~SDCustomCollapse() {}
  double calcSD(const arma::vec& f,
                const arma::vec& t) {
    mat ftMat = CollapseTree(this->edgeMat, this->tipNum, f, t);
    vec fnew = ftMat.col(2);
    vec tnew = ftMat.col(3);
    return func(fnew, tnew);
  }
};


#endif
