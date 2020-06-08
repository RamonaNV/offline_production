#include "simclasses/I3CorsikaWeight.h"

static const unsigned i3corsika_weight_version_ = 0;

std::ostream& I3CorsikaWeight::Print(std::ostream& oss) const{
  I3Position pos = primary.GetPos();
  oss << "[ I3CorsikaWeight ["
      << primary.GetMinorID() << " " << primary.GetTypeString() << " "
      << "(" << primary.GetZenith()/I3Units::degree << "deg, "
      << primary.GetAzimuth()/I3Units::degree << "deg) "
      << primary.GetEnergy()/I3Units::GeV << "GeV], "
      << bias
      << ", weight: " << std::fixed << std::setprecision(4) << weight
      << ", max_x: "  << std::fixed << std::setprecision(4) << max_x << " ]";
  return oss;
}

template <class Archive> 
void I3CorsikaWeight::serialize(Archive& ar, unsigned version){
  if (version>i3corsika_weight_version_){
    log_fatal("Attempting to read version %u from file but running version %u of I3Position class.",
	      version,i3position_version_);
  }
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("primary", primary);
  ar & make_nvp("bias",bias);
  ar & make_nvp("weight", weight);
  ar & make_nvp("max_x", max_x);

}

std::ostream& operator<<(std::ostream& oss, const I3CorsikaWeight& w){
  return(w.Print(oss));
}

I3_CLASS_VERSION(I3CorsikaWeight, i3corsika_weight_version_);
I3_SERIALIZABLE(I3CorsikaWeight);
