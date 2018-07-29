#include "utilities.h"

bool isEqualStr(std::string& str1, std::string str2) {
  return str1.compare(str2) == 0;
}

double similarityToDistance(const double distance) {
  return 1.0 - std::abs(distance);
}

