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
  else if (isEqualStr(sdName, "SimJaccard")) {
    sdfunc = std::make_shared<SimJaccard>();
  }
  else if (isEqualStr(sdName, "DistHamming")) {
    sdfunc = std::make_shared<DistHamming>();
  }
  else if (isEqualStr(sdName, "DistManhattan")) {
    sdfunc = std::make_shared<DistManhattan>();
  }
  else if (isEqualStr(sdName, "DistEuclidean")) {
    sdfunc = std::make_shared<DistEuclidean>();
  }
  else if (isEqualStr(sdName, "DistMinkowski")) {
    int p = 2;
    if (arguments.containsElementNamed("p")) {
      p = Rcpp::as<int>(arguments["p"]);
    } else {}
    sdfunc = std::make_shared<DistMinkowski>(p);
  }
  else if (isEqualStr(sdName, "SimMIBin")) {
    sdfunc = std::make_shared<SimMIBin>();
  }
  else if (isEqualStr(sdName, "SimMIConti")) {
    uword bin = 10;
    if (arguments.containsElementNamed("bin")) {
      bin = Rcpp::as<uword>(arguments["bin"]);
    } else {}
    sdfunc = std::make_shared<SimMIConti>(bin);
  }
  else if (isEqualStr(sdName, "custom")) {
    SEXP func_ = arguments["func"];
    funcPtr func = *Rcpp::XPtr<funcPtr>(func_);
    sdfunc = std::make_shared<SDCustom>(func);
  } else {}

  return sdfunc;
}
