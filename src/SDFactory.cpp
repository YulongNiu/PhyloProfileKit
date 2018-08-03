#include "SDFactory.h"
#include "SD.h"
#include "utilities.h"

std::shared_ptr<SDmeasure> SDFactory::createSDFunc(Rcpp::List &attrs,
                                                   Rcpp::List &arguments) {

  std::string sdName = attrs["method"];
  std::shared_ptr<SDmeasure> sdfunc = NULL;

  if (isEqualStr(sdName, "SimCor")) {
    sdfunc = std::make_shared<SimCor>();
  }
  else if (isEqualStr(sdName, "SimCorCollapse")) {
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<SimCorCollapse>(edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "SimJaccard")) {
    sdfunc = std::make_shared<SimJaccard>();
  }
  else if (isEqualStr(sdName, "SimJaccardCollapse")) {
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<SimJaccardCollapse>(edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "DistHamming")) {
    sdfunc = std::make_shared<DistHamming>();
  }
  else if (isEqualStr(sdName, "DistHammingCollapse")) {
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<DistHammingCollapse>(edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "DistManhattan")) {
    sdfunc = std::make_shared<DistManhattan>();
  }
  else if (isEqualStr(sdName, "DistManhattanCollapse")) {
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<DistManhattanCollapse>(edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "DistEuclidean")) {
    sdfunc = std::make_shared<DistEuclidean>();
  }
  else if (isEqualStr(sdName, "DistEuclideanCollapse")) {
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<DistEuclideanCollapse>(edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "DistMinkowski")) {
    uword p = 3;
    if (arguments.containsElementNamed("p")) {
      p = Rcpp::as<uword>(arguments["p"]);
    } else {}
    sdfunc = std::make_shared<DistMinkowski>(p);
  }
  else if (isEqualStr(sdName, "DistMinkowskiCollapse")) {
    uword p = 3;
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    if (arguments.containsElementNamed("p")) {
      p = Rcpp::as<uword>(arguments["p"]);
    } else {}
    sdfunc = std::make_shared<DistMinkowskiCollapse>(p, edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "SimMIBin")) {
    sdfunc = std::make_shared<SimMIBin>();
  }
  else if (isEqualStr(sdName, "SimMIBinCollapse")) {
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<SimMIBinCollapse>(edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "SimMIConti")) {
    uword bin = 10;
    if (arguments.containsElementNamed("bin")) {
      bin = Rcpp::as<uword>(arguments["bin"]);
    } else {}
    sdfunc = std::make_shared<SimMIConti>(bin);
  }
  else if (isEqualStr(sdName, "SimMIContiCollapse")) {
    uword bin = 10;
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    if (arguments.containsElementNamed("bin")) {
      bin = Rcpp::as<uword>(arguments["bin"]);
    } else {}
    sdfunc = std::make_shared<SimMIContiCollapse>(bin, edgeMat, tipNum);
  }
  else if (isEqualStr(sdName, "SDCustom")) {
    SEXP func_=arguments["func"];
    funcPtr func = *Rcpp::XPtr<funcPtr>(func_);
    sdfunc = std::make_shared<SDCustom>(func);
  }
  else if (isEqualStr(sdName, "SDCustomCollapse")) {
    SEXP func_=arguments["func"];
    funcPtr func = *Rcpp::XPtr<funcPtr>(func_);
    mat edgeMat = arguments["edgeMat"];
    uword tipNum = arguments["tipNum"];
    sdfunc = std::make_shared<SDCustomCollapse>(func, edgeMat, tipNum);
  }
  else {}

  return sdfunc;
}
