#include "SDFactory.h"
#include "SD.h"
#include "Utility.h"

std::shared_ptr<SDmeasure> SDFactory::createSDFunc(Rcpp::List &attrs,
                                                   Rcpp::List &arguments) {

  using namespace utility;
  std::string sdName = attrs["method"];
  std::shared_ptr<SDmeasure> sdfunc = NULL;

  if (isEqualStr(sdName, "distsum")) {
    sdfunc = std::make_shared<DistSum>();
  } else if (isEqualStr(sdName, "simjaccard")) {
    sdfunc = std::make_shared<SimJaccard>();
  } else if (isEqualStr(sdName, "distminkowski")) {
    int p = 2;
    if (arguments.containsElementNamed("p")) {
      p = Rcpp::as<int>(arguments["p"]);
    }
    sdfunc = std::make_shared<DistMinkowski>(p);
  }
  else {}

  return sdfunc;
}
