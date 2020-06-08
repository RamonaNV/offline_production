#ifndef I3CORSIKA_INFO_H_INCLUDED
#define I3CORSIKA_INFO_H_INCLUDED

#include "icetray/I3DefaultName.h"
#include "dataclasses/physics/I3Particle.h"

static const unsigned i3corsika_info_version_ = 0;

struct I3CorsikaInfo : public I3FrameObject{

  uint32_t run_id;
  uint32_t n_events;
  I3Particle::ParticleType primary_type;
  uint32_t atmosphere;
  uint32_t oversampling;
  double cylinder_height;
  double cylinder_radius;
  double min_zenith;
  double max_zenith;
  double min_energy;
  double max_energy;
  double power_law_index;

  I3CorsikaInfo();

  std::ostream& Print(std::ostream&) const override;

  friend class icecube::serialization::access;

  template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

std::ostream& operator<<(std::ostream& oss, const I3CorsikaInfo& n);

I3_POINTER_TYPEDEFS(I3CorsikaInfo);
I3_DEFAULT_NAME(I3CorsikaInfo);
I3_CLASS_VERSION(I3CorsikaInfo, i3corsika_info_version_);

#endif //I3CORSIKA_INFO_H_INCLUDED
