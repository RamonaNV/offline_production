#include <icetray/serialization.h>
#include <simclasses/I3CorsikaInfo.h>

I3CorsikaInfo::I3CorsikaInfo() :
  I3FrameObject(), run_id(0), n_events(0), primary_type(I3Particle::unknown),
  atmosphere(0), oversampling(0), cylinder_height(NAN), cylinder_radius(NAN), 
  min_zenith(NAN), max_zenith(NAN), min_energy(NAN), max_energy(NAN), 
  power_law_index(NAN) {}

template <class Archive>
void 
I3CorsikaInfo::serialize(Archive& ar, unsigned version)
{
  if (version>i3corsika_info_version_){
    log_fatal("Attempting to read version %u from file but running version %u of I3CorsikaInfo class.",version,i3corsika_info_version_);
  }
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("run_id", run_id);
  ar & make_nvp("n_events",n_events);  
  ar & make_nvp("primary_type",primary_type);
  ar & make_nvp("atmosphere",atmosphere);  
  ar & make_nvp("oversampling",oversampling);  
  ar & make_nvp("cylinder_height",cylinder_height);
  ar & make_nvp("cylinder_radius",cylinder_radius);
  ar & make_nvp("min_zenith",min_zenith);
  ar & make_nvp("max_zenith",max_zenith);  
  ar & make_nvp("min_energy",min_energy);  
  ar & make_nvp("max_energy",max_energy);
  ar & make_nvp("power_law_index",power_law_index);
}

std::ostream& I3CorsikaInfo::Print(std::ostream& oss) const{  
  oss << "[ I3CorsikaInfo" << std::endl
      << "            run_id : " << run_id << std::endl
      << "          n_events : " << n_events << std::endl
      << "      primary_type : " << i3particle_type_string(primary_type) << std::endl
      << "        atmosphere : " << atmosphere << std::endl    
      << "      oversampling : " << oversampling << std::endl
      << "   cylinder_height : " << cylinder_height << std::endl
      << "   cylinder_radius : " << cylinder_radius << std::endl
      << "        min_zenith : " << min_zenith << std::endl
      << "        max_zenith : " << max_zenith << std::endl
      << "        min_energy : " << min_energy << std::endl    
      << "        max_energy : " << max_energy << std::endl
      << "   power_law_index : " << power_law_index << std::endl
      << "]" ;
  return oss;
}

std::ostream& operator<<(std::ostream& oss, const I3CorsikaInfo& n){
  return(n.Print(oss));
}

I3_SERIALIZABLE(I3CorsikaInfo);
