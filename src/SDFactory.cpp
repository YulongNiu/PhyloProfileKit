#include "SDFactory.h"
#include "SD.h"
#include "utilities.h"

std::shared_ptr<SDmeasure> SDFactory::createSDFunc(Rcpp::List &attrs,
                                                   Rcpp::List &arguments) {

  std::string sdName = attrs["method"];
  std::shared_ptr<SDmeasure> sdfunc = NULL;

  if (isEqualStr(sdName, "SimJaccard")) {
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
  } else {}

  return sdfunc;
}
