#ifndef I3SHOWERBIAS_H_INCLUDED
#define I3SHOWERBIAS_H_INCLUDED

#include "icetray/I3DefaultName.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/physics/I3Particle.h"

struct I3ShowerBias : public I3FrameObject{

  enum BiasParticleType {
    Mu = 0,
    NuMu,
    NuE
  };

  BiasParticleType type;
  double target;  

  I3ShowerBias(){ type=BiasParticleType(-1);target=NAN;}
  
  I3ShowerBias(BiasParticleType t1, double t2){
    type=t1; target=t2;
  }

  static std::string bias_particle_type_string(BiasParticleType);

  std::ostream& Print(std::ostream&) const override;

  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive& ar, unsigned version);
};

std::ostream& operator<<(std::ostream& oss, const I3ShowerBias& d);

typedef I3Map<I3ParticleID, I3ShowerBias> I3ShowerBiasMap;

I3_DEFAULT_NAME(I3ShowerBiasMap);
I3_POINTER_TYPEDEFS(I3ShowerBiasMap);

#endif
