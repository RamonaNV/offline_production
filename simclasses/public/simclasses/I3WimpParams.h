/**
 * copyright  (C) 2012
 * the icecube collaboration
 * @version $Id: I3WimpParams.h xxxxx 2012-06-01 17:01:33Z mzoll $
 * @file I3WimpParams.h
 * @date $Date: 2012-12-20$
 * @author mzoll <marcel.zoll@fysik.su.se>
 * This is a archive class that saves important information about WIMP events that are read by the WimpSim-Reader
 */

#ifndef I3WIMPPARAMS_H_INCLUDED
#define I3WIMPPARAMS_H_INCLUDED

#define WIMPSIMPARAMS_H_WimpSim_DecayChannel                                  \
    (unknown)(down)(up)(strange)(charm)(bottom)(top)(gluon)(w)(z)(muon)(tau)      \
    (nue)(numu)(nutau)(KKDM)
    
#define WIMPSIMPARAMS_H_WimpSim_PythiaParticleTypes                           \
    (Down(DownBar)(Up)(UpBar)(Strange)(StrangeBar)(Charm)(CharmBar)(Bottom)       \
    (BottomBar)(Tau)(TauBar)(Gluon)(WPlus)(WMinus)(ZNull)(Znull)(MuMinus)(MuPlus) \
    (TauMinus)(TauPlus)(NuE)(NuEBar)(NuMu)(NuMuBar)(NuTau)(NuTauBar)

#define I3WIMPPARAMS_H_I3WimpParams_Source                                        \
    (UNKNOWN)(SUN)(EARTH)
static const unsigned i3wimpparams_version_ = 0;

#include <string>
//#include "wimpsim-reader/WimpSimTools.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3DefaultName.h"
#include "icetray/serialization.h"
#include "dataclasses/Utility.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/I3Time.h"


/// declaration of the particle types that are included in the pythia code of WimpSim
namespace PythiaParticleTypes{
  /// enum of the Paricle Types used by pythia
  enum ParticleTypes {
    Down = 1,
    DownBar = -1,
    Up = 2,
    UpBar = -2,
    Strange = 3,
    StrangeBar = -3,
    Charm = 4,
    CharmBar = -4,
    Bottom = 5,
    BottomBar = -5,
    Tau = 6,
    TauBar = -6,
    Gluon = 21,
    WPlus = 24,
    WMinus = -24,
    ZNull = 23,
    Znull = -23,
    MuMinus = 13,
    MuPlus = -13,
    TauMinus = 15,
    TauPlus = -15,
    NuE = 12,
    NuEBar = -12,
    NuMu = 14,
    NuMuBar = -14,
    NuTau = 16,
    NuTauBar = -16
  };
};

/// Declaration of WimpSim internal conventions
namespace WimpSim{
  /// enum of the Decay Channels used by WimpSim
  enum DecayChannel {
    unknown = 0,
    down = 1,
    up = 2,
    strange = 3,
    charm = 4,
    bottom = 5,
    top = 6,
    gluon = 7,
    w = 8,
    z = 9,
    muon = 10,
    tau = 11,
    nue = 12,
    numu = 13,
    nutau = 14,
    KKDM = 100
  };
  /// enum of the Particle Types used by WimpSim
  enum ParticleType{
    error = 0,
    e_plus = 10,
    e_minus = 11,
    nu_e = 15,
    nu_e_bar = 16,
    mu_plus = 20,
    mu_minus = 21,
    nu_mu = 25,
    nu_mu_bar = 26,
    tau_plus = 30,
    tau_minus = 31,
    nu_tau = 35,
    nu_tau_bar = 36,
    hshow = 100,
  };
  
  
    /// enumerator of the Source the Wimp is comming from
  enum SourceType {
    UNKNOWN = 0,
    SUN = 1,
    EARTH = 2
  };
  
  
  /** @brief converts WimpSim-Particle string to Particle
    @param particletype : particle from wimpsim as string
    @return ParticleType as enum
  */
  ParticleType StringToParticleType (std::string particletype);
  
  /** @brief converts WimpSim-Particle string to I3Particle conventions
    @param particletype : particle from wimpsim as string
    @return the correct I3ParticleType
    */
  I3Particle::ParticleType WimpSimToI3ParticleType (std::string particletype);
  
  /** @brief converts WimpSim-Particle to I3Particle conventions
    @param particletype : ParticleType from WimpSim
    @return the correct I3ParticleType
    */
  I3Particle::ParticleType ParticleTypeToI3ParticleType (ParticleType particletype);
  
  /** @brief converts I3Particle string to WimpSim-Particle conventions
    @param particletype : ParticleType from I3Particle
    @return the correct WimpSimParticleType
    */
  ParticleType I3ParticleTypeToParticleType (I3Particle::ParticleType particletype);
  
  /** @brief is this particle a lepton(any neutrino, or charged lepton)
    @param particle : particle from wimpsim
    @return it is a lepton
  */  
  bool IsLepton(ParticleType particle);
  
  /** @brief is this particle a neutrino
    @param particle : particle from wimpsim
    @return it is a neutrino
  */  
  bool IsNeutrino(ParticleType particle);
  
  /** @brief is this particle a charged lepton
    @param particle : particle from wimpsim
    @return it is a charged lepton
    */
  bool IsChargedLepton(ParticleType particle);
  
  /** @brief convert Decay Channel int type to enum
    @param channelnbr : a plain int
    @return DecayChannel as enum type
    */
  DecayChannel ConvertToDecayChannel(int channelnbr);
  
  /** @brief Convert enum type Decay Channel to plain int
    @param channel : from WimpSim
    @return DecayChannel as an int
  */
  unsigned int DecayChannelToInt (DecayChannel channel);
  
   /** @brief Convert Source as a string to a int
   * @param sourcestr : from WimpSim
   * @return 1 for Sun
   *       2 for Earth
   *       0 for default
   */
  unsigned int SourceStringToInt(std::string sourcestr);
}


//=======================================================

/**
 * The basic WimpParams class for WimpSimReader for IceCube.
 */
class I3WimpParams : public I3FrameObject {
private: // these are the saved parameters
  /// The origin of the WIMP
  WimpSim::SourceType source_;
  /// The mass of the WIMP
  double mass_;
  /// The simulated decay Channel of the WIMP
  WimpSim::DecayChannel channel_;
  /// The weight of the neutrino
  double nu_weight_;
  /// The weight of the inice lepton
  double lep_weight_;
  /// The weight of the inice hadron
  double had_weight_;
  /// The generated volume for this event
  double vgen_;
  /// The neutrino generated/effective area for this event
  double aproj_;
  /// The time for this event
  I3Time time_;

public:
  /**
  * stdConstructor; set everything to NAN or equivalent
  */
  I3WimpParams();

  /**
  * stdDestructor
  */
  virtual ~I3WimpParams();

  /**
  * Constructor with enum
  */
  I3WimpParams(WimpSim::SourceType source,
               double mass,
               WimpSim::DecayChannel channel,
               double nu_weight,
               double lep_weight,
               double had_weight,
               double vgen,
               double aproj,
               I3Time time);

  /**
  * Constructor with basic data-types only
  */
  I3WimpParams(std::string source,
               double mass,
               int channel,
               double nu_weight,
               double lep_weight,
               double had_weight,
               double vgen,
               double aproj,
               I3Time time);
  
  std::ostream& Print(std::ostream&) const override;

  /////////// Getters /////////
  /**
  * Get the source of the Wimp as enum representation
  */
  WimpSim::SourceType GetSource() const;
  /**
  * Get the source of the Wimp as string representation
  */
  std::string GetSourceString() const;
  /**
  * Get the mass
  */
  double GetMass() const;
  /**
  * Get the channel
  */
  WimpSim::DecayChannel GetChannel() const;
  /**
  * Get the channel int representation
  */
  int GetChannelInt() const;
  /**
  * Get the primary neutrino weight as natively saved by WIMPsim
  */
  double GetNuWeight() const;
  /**
  * Get the leptonic weight as natively saved by WIMPsim
  */
  double GetLepWeight() const;
  /**
  * Get the hadronic weight as natively saved by WIMPsim
  */
  double GetHadWeight() const;
  /**
  * Get the leptonic weight multiplied by E47
  */
  double GetLepWeightE47() const;
  /**
  * Get the hadronic weight multiplied by E47
  */
  double GetHadWeightE47() const;
  /**
  * Get the generated Volume
  */
  double GetVgen() const;
  /**
  * Get the generated Volume
  */
  double GetAproj() const;
  /**
  * Get the Time
  */
  I3Time GetTime() const;

  /////////// Setters ///////
  /**
  * Set the source as SourceType
  */
  void SetSource(const WimpSim::SourceType source);
  /**
  * Set the source as String
  */
  void SetSource(const std::string source);
  /**
  * Set the source as a uint
  */
  void SetSource(const uint source);
  /**
  * Set the mass
  */
  void SetMass(const double mass);
  /**
  * Set the channel
  */
  void SetChannel(const WimpSim::DecayChannel channel);
  /**
  * Set the primary neutrino weight as natively saved by WIMPsim
  */
  void SetNuWeight(const double nu_weight);
  /**
  * Set the leptonic weight as natively saved by WIMPsim
  */
  void SetLepWeight(const double lep_weight);
  /**
  * Set the hadronic weight as natively saved by WIMPsim
  */
  void SetHadWeight(const double had_weight);
  /**
  * Set the leptonic weight already multiplied by E47
  */
  void SetLepWeightE47(const double lep_weight);
  /**
  * Set the hadronic weight already multiplied by E47
  */
  void SetHadWeightE47(const double had_weight);
  /**
  * Set the generated Volume
  */
  void SetVgen(const double vgen);
  /**
  * Set the generated Volume
  */
  void SetAproj(const double aproj);
  /**
  * Set the Time from a I3Time
  */
  void SetTime(const I3Time time);
  /**
    * Set the time as MJD-double
    */
  void SetTime(const double mjd);

private: // This is I3Tray python and serialization stuff
  /**
  * Serialization to I3
  */
  friend class icecube::serialization::access;
  /**
  * make Archive-type in I3Tray
  */
  template <class Archive>
  /**
  * serialization away
  */
  void serialize(Archive& ar, unsigned version);
};

/**
* Dump the object into a string-stream (this is bound to python as string-cast)
*/
std::ostream& operator<<(std::ostream& oss, const I3WimpParams& d);

I3_POINTER_TYPEDEFS(I3WimpParams);
I3_CLASS_VERSION(I3WimpParams, i3wimpparams_version_);
I3_DEFAULT_NAME(I3WimpParams);

#endif //I3WIMPPARAMS_H_INCLUDED
