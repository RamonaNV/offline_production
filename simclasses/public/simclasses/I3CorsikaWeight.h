#include "simclasses/I3ShowerBias.h"
#include "dataclasses/physics/I3Particle.h"

struct I3CorsikaWeight : public I3FrameObject
{
  I3Particle primary;
  I3ShowerBias bias;
  double weight;
  double max_x;

  std::ostream& Print(std::ostream&) const override;

private:
  friend class icecube::serialization::access;
  
  template <class Archive>
  void serialize(Archive& ar, unsigned version);
};

std::ostream& operator<<(std::ostream& oss, const I3CorsikaWeight& d);

I3_POINTER_TYPEDEFS(I3CorsikaWeight);
