/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef I3MCPE_H_INCLUDED
#define I3MCPE_H_INCLUDED

#include <vector>
#include <icetray/I3Logging.h>
#include <icetray/serialization.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/physics/I3ParticleID.h>
#include <ostream>

static const unsigned i3mcpe_version_ = 1;

/**
 * @brief I3MCPE struct that stores the photon arrival time 
 * (i.e.PE creation time), number of PE (for binning), and
 * the IDs of the particle that created this.
 */

struct I3MCPE {

  /**
   * ID of the I3Particle that created this PE
   */
  I3ParticleID ID;

  /**
   * Creation time of PE (photon arrival time)
   */ 
  double time;

  /**
   * Number of PEs this object represents.
   * Used for binning.
   */
  uint32_t npe;

  SET_LOGGER("I3MCPE");

  bool operator==(const I3MCPE& rhs) const {
    return time == rhs.time
    && npe == rhs.npe
    && ID == rhs.ID;
  }
	
  // default constructor for noise generators
  // I3MCPE():npe(0)
  // {
  //   ID.majorID = 0;
  //   ID.minorID = 0;
  // }
  
  I3MCPE(const uint32_t npe_ = 0, const double time_ = 0):
  time(time_), npe(npe_)
  {
    ID.majorID = 0;
    ID.minorID = 0;
  }

  // constructor for hit makers
  // this just sets the major and minor IDs accordingly
  I3MCPE(const I3ParticleID& p, const uint32_t n_pe = 0, const double pe_time = 0):
  ID(p),time(pe_time),npe(n_pe){}

  I3MCPE(const uint64_t major_ID, const int32_t minor_ID, const uint32_t n_pe = 0, const double pe_time = 0):
  time(pe_time), npe(n_pe)
  {
    ID.majorID = major_ID;
    ID.minorID = minor_ID;
  }
  
  operator I3ParticleID() const{ return ID; }
  
  std::ostream& Print(std::ostream&) const;
  
private:
  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive & ar, const unsigned version)
  {
    if (version>i3mcpe_version_)
      log_fatal("Attempting to read version %u from file but running version %u of I3MCPE class.",
                version,i3mcpe_version_);
    if(version == 0){
      float t(0.);
      ar & make_nvp("time",t);
      time = t;
    }else{
      ar & make_nvp("time",time);
    }
    ar & make_nvp("npe",npe);
    ar & make_nvp("major_ID",ID.majorID);
    ar & make_nvp("minor_ID",ID.minorID);
  }

};

I3_CLASS_VERSION(I3MCPE,i3mcpe_version_);

typedef std::vector<I3MCPE> I3MCPESeries;
typedef I3Map<OMKey, I3MCPESeries > I3MCPESeriesMap;

std::ostream& operator<<(std::ostream&, const I3MCPE&);
std::ostream& operator<<(std::ostream&, const I3MCPESeries&);

I3_POINTER_TYPEDEFS(I3MCPE);
I3_POINTER_TYPEDEFS(I3MCPESeries);
I3_POINTER_TYPEDEFS(I3MCPESeriesMap);

#endif
