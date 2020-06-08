#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include "MuonGun/WeightCalculator.h"

namespace I3MuonGun {

std::string get_tabledir();

BundleModel load_model(const std::string &base);

}

#endif  // COMMON_H_INCLUDED
