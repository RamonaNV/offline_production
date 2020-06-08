#include <icetray/serialization.h>
#include <simclasses/I3ShowerBias.h>


std::string I3ShowerBias::bias_particle_type_string(BiasParticleType t)
{
  switch (t) {
  case Mu:   return "Mu";  break;
  case NuMu: return "NuMu";break;        
  case NuE:  return "NuE"; break;    
  default:   return "Unknown";
  }
}

std::ostream& I3ShowerBias::Print(std::ostream& oss) const{  
  oss << "I3ShowerBias("
      << bias_particle_type_string(type) << ", "
      << std::fixed << std::setprecision(6) << target << ")";
  return oss;
}

std::ostream& operator<<(std::ostream& oss, const I3ShowerBias& n){
  return(n.Print(oss));
}

template <class Archive>
void 
I3ShowerBias::serialize(Archive& ar, unsigned version)
{
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("type",type);
  ar & make_nvp("target",target);
}

I3_SERIALIZABLE(I3ShowerBias);
I3_SERIALIZABLE(I3ShowerBiasMap);
