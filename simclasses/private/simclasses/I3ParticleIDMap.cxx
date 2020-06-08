#include <simclasses/I3ParticleIDMap.hpp>
#include <ostream>
#include "I3SequenceOpOStream.h"

I3_SERIALIZABLE(I3ParticleIDMap);

std::ostream& operator<<(std::ostream& os, const ParticlePulseIndexMap& ppim){
  for(const auto& particlePair : ppim){
    os << "  " << particlePair.first << ": [";
    bool first=true;
    for(auto idx : particlePair.second){
      if(first)
        first=false;
      else
        os << ", ";
      os << idx;
    }
    os << "]\n";
  }
  return os;
}
