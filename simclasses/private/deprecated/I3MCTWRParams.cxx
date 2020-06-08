#include <icetray/serialization.h>
#include "icetray/OMKey.h"

static const unsigned i3mctwrparams_version_ = 1;

struct I3MCTWRParams{

  int stop_delay;
  int DMADD_thresh;
  int TWR_thresh;
  double rel_sens;
  int wf_type;
  double afterpulse_prob;
  double afterpulse_time;
  double noise_rate;
  double amplitude;
  double cable_delay;
  bool optical;
  double peArea;

  virtual ~I3MCTWRParams(){};
  
  std::ostream& Print(std::ostream&) const;

  template <class Archive> void serialize(Archive & ar, unsigned version);
};

typedef std::map<OMKey, I3MCTWRParams> I3MCTWRParamsMap;


I3_CLASS_VERSION(I3MCTWRParams, i3mctwrparams_version_);

template <class Archive>
  void I3MCTWRParams::serialize(Archive& ar, unsigned version)
  {
    ar & make_nvp("stop_delay",stop_delay);
    ar & make_nvp("DMADD_thresh",DMADD_thresh);
    ar & make_nvp("TWR_thresh",TWR_thresh);
    ar & make_nvp("rel_sens",rel_sens);
    ar & make_nvp("wf_type",wf_type);
    ar & make_nvp("afterpulse_prob",afterpulse_prob);
    ar & make_nvp("afterpulse_t",afterpulse_time);
    ar & make_nvp("noise_rate",noise_rate);
    ar & make_nvp("amplitude",amplitude);
    ar & make_nvp("cable_delay",cable_delay);
    ar & make_nvp("optical",optical);
    if(version > 0){
      ar & make_nvp("peArea",peArea);    
    }
  }

I3_SERIALIZABLE(I3MCTWRParams);

std::ostream& I3MCTWRParams::Print(std::ostream& os) const{
  os << "[I3MCTWRParams:\n"
     << "        Stop Delay: " << stop_delay << '\n'
     << "   DMADD Threshold: " << DMADD_thresh << '\n'
     << "     TWR Threshold: " << TWR_thresh << '\n'
     << "    Relative Sens.: " << rel_sens << '\n'
     << "     Waveform Type: " << wf_type << '\n'
     << "  Afterpulse Prob.: " << afterpulse_prob << '\n'
     << "   Afterpulse Time: " << afterpulse_time << '\n'
     << "        Noise Rate: " << noise_rate << '\n'
     << "         Amplitude: " << amplitude << '\n'
     << "       Cable Delay: " << cable_delay << '\n'
     << "           Optical: " << optical << '\n'
     << "         P.E. Area: " << peArea << "\n]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const I3MCTWRParams& p){
  return(p.Print(os));
}

