/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef I3MCPULSE_H_INCLUDED
#define I3MCPULSE_H_INCLUDED

#include <algorithm>
#include <iterator>
#include <vector>
#include <boost/foreach.hpp>
#include <icetray/I3Logging.h>
#include <icetray/serialization.h>
#include <dataclasses/I3Map.h>

static const unsigned i3mcpulse_version_ = 1;

/**
 * @brief I3MCPulse struct that stores the time, charge,
 * and source of a PMT pulse.
 */

#define I3MCPULSE_H_I3MCPulse_PulseSource      \
  (UNKNOWN)(PE)(RANDOM)(AFTER_PULSE)(PRE_PULSE)\
  (ELASTIC_LATE_PULSE)(INELASTIC_LATE_PULSE)   \
  (EARLY_AFTER_PULSE)

struct I3MCPulse {

  enum PulseSource{
    UNKNOWN = 0,
    PE = 10,
    RANDOM = 20,
    AFTER_PULSE = 30,
    PRE_PULSE = 40,
    ELASTIC_LATE_PULSE = 50,
    INELASTIC_LATE_PULSE = 60,
    EARLY_AFTER_PULSE = 70,
    CROSSTALK_PULSE = 80
  };

  double time;
  float charge;
  PulseSource source;

  SET_LOGGER("I3MCPulse");

  bool operator==(const I3MCPulse& rhs) const {
    return time == rhs.time
    && charge == rhs.charge
    && source == rhs.source;
  }

  I3MCPulse(float t, float c = 1.0 , PulseSource s = PE) :
    time(t), charge(c), source(s) {};

  I3MCPulse(){};
  
  std::ostream& Print(std::ostream&) const;
  
private:
  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive & ar, const unsigned version)
  {
    if (version>i3mcpulse_version_)
      log_fatal("Attempting to read version %u from file but running version %u of I3MCPulse class.",
                version,i3mcpulse_version_);
    if(version == 0){
      float t(0.);
      ar & make_nvp("time",t);
      time = t;
    }else{
      ar & make_nvp("time",time);
    }
    ar & make_nvp("charge",charge);
    ar & make_nvp("source",source);
  }

};

I3_CLASS_VERSION(I3MCPulse,i3mcpulse_version_);

typedef std::vector<I3MCPulse> I3MCPulseSeries;
typedef I3Map<OMKey, I3MCPulseSeries > I3MCPulseSeriesMap;

std::ostream& operator<<(std::ostream&, const I3MCPulse&);
std::ostream& operator<<(std::ostream&, const I3MCPulseSeries&);

I3_POINTER_TYPEDEFS(I3MCPulseSeries);
I3_POINTER_TYPEDEFS(I3MCPulseSeriesMap);

#endif
