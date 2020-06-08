/**
 * copyright  (C) 2012
 * the icecube collaboration
 * @version $Id: $
 * @file I3WimpParams.cxx
 * @date $Date: 2012-11-18$
 * @author mzoll <marcel.zoll@fysik.su.se>
 * This is a archive class that saves important information about WIMP Events that are read by the WimpSim-Reader
 */

#include <cmath>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "icetray/serialization.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"
#include "simclasses/I3WimpParams.h"



//_______________________________________________________________________
WimpSim::DecayChannel WimpSim::ConvertToDecayChannel (int channelnbr) {
  if (channelnbr == 0)
    return unknown;
  else if (channelnbr == 1)
    return down;
  else if (channelnbr == 2)
    return up;
  else if (channelnbr == 3)
    return strange;
  else if (channelnbr == 4)
    return charm;
  else if (channelnbr == 5)
    return bottom;
  else if (channelnbr == 6)
    return top;
  else if (channelnbr == 7)
    return gluon;
  else if (channelnbr == 8)
    return w;
  else if (channelnbr == 9)
    return z;
  else if (channelnbr == 10)
    return muon;
  else if (channelnbr == 11)
    return tau;
  else if (channelnbr == 12)
    return nue;
  else if (channelnbr == 13)
    return numu;
  else if (channelnbr == 14)
    return nutau;
  else if (channelnbr == 100)
    return KKDM;
  else
    return unknown;
}
//_______________________________________________________________________
unsigned int WimpSim::DecayChannelToInt (WimpSim::DecayChannel channel) {
  return channel;
}
//_______________________________________________________________________
I3Particle::ParticleType WimpSim::WimpSimToI3ParticleType (std::string particletype) { //NOTE this is not a complete list
  if (particletype == "e+")
    return I3Particle::EPlus;
  else if (particletype == "e-")
    return I3Particle::EMinus;
  else if (particletype == "mu+")
    return I3Particle::MuPlus;
  else if (particletype == "mu-")
    return I3Particle::MuMinus;
  else if (particletype == "tau+")
    return I3Particle::TauPlus;
  else if (particletype == "tau-")
    return I3Particle::TauMinus;
  else if (particletype == "nu_e")
    return I3Particle::NuE;
  else if (particletype == "nu_e~")
    return I3Particle::NuEBar;
  else if (particletype == "nu_mu")
    return I3Particle::NuMu;
  else if (particletype == "nu_mu~")
    return I3Particle::NuMuBar;
  else if (particletype == "nu_tau")
    return I3Particle::NuTau;
  else if (particletype == "nu_tau~")
    return I3Particle::NuTauBar;
  else if (particletype == "hshow")
    return I3Particle::Hadrons;
  else
    return I3Particle::unknown;
}
//_______________________________________________________________________
WimpSim::ParticleType WimpSim::StringToParticleType (std::string particletype) {
  if (particletype == "e+")
    return e_plus;
  else if (particletype == "e-")
    return e_minus;
  else if (particletype == "mu+")
    return mu_plus;
  else if (particletype == "mu-")
    return mu_minus;
  else if (particletype == "tau+")
    return tau_plus;
  else if (particletype == "tau-")
    return tau_minus;
  else if (particletype == "nu_e")
    return nu_e;
  else if (particletype == "nu_e~")
    return nu_e_bar;
  else if (particletype == "nu_mu")
    return nu_mu;
  else if (particletype == "nu_mu~")
    return nu_mu_bar;
  else if (particletype == "nu_tau")
    return nu_tau;
  else if (particletype == "nu_tau~")
    return nu_tau_bar;
  else if (particletype == "hshow")
    return hshow;
  else
    return error;
}
//_______________________________________________________________________
I3Particle::ParticleType WimpSim::ParticleTypeToI3ParticleType (ParticleType particletype) { //NOTE this is not a complete list
  if (particletype == e_plus)
    return I3Particle::EPlus;
  else if (particletype == e_minus)
    return I3Particle::EMinus;
  else if (particletype == mu_plus)
    return I3Particle::MuPlus;
  else if (particletype == mu_minus)
    return I3Particle::MuMinus;
  else if (particletype == tau_plus)
    return I3Particle::TauPlus;
  else if (particletype == tau_minus)
    return I3Particle::TauMinus;
  else if (particletype == nu_e)
    return I3Particle::NuE;
  else if (particletype == nu_e_bar)
    return I3Particle::NuEBar;
  else if (particletype == nu_mu)
    return I3Particle::NuMu;
  else if (particletype == nu_mu_bar)
    return I3Particle::NuMuBar;
  else if (particletype == nu_tau)
    return I3Particle::NuTau;
  else if (particletype == nu_tau_bar)
    return I3Particle::NuTauBar;
  else if (particletype == hshow)
    return I3Particle::Hadrons;
  else
    return I3Particle::unknown;
}
//_______________________________________________________________________
WimpSim::ParticleType WimpSim::I3ParticleTypeToParticleType (I3Particle::ParticleType particletype) { //NOTE this is not a complete list
  if (particletype == I3Particle::EPlus)
    return e_plus;
  else if (particletype == I3Particle::EMinus)
    return e_minus;
  else if (particletype == I3Particle::MuPlus)
    return mu_plus;
  else if (particletype == I3Particle::MuMinus)
    return mu_minus;
  else if (particletype == I3Particle::TauPlus)
    return tau_plus;
  else if (particletype == I3Particle::TauMinus)
    return tau_minus;
  else if (particletype == I3Particle::NuE)
    return nu_e;
  else if (particletype == I3Particle::NuEBar)
    return nu_e_bar;
  else if (particletype == I3Particle::NuMu)
    return nu_mu;
  else if (particletype == I3Particle::NuMuBar)
    return nu_mu_bar;
  else if (particletype == I3Particle::NuTau)
    return nu_tau;
  else if (particletype == I3Particle::NuTauBar)
    return nu_tau_bar;
  else if (particletype == I3Particle::Hadrons)
    return hshow;
  else
    return error;
}
//_______________________________________________________________________
bool WimpSim::IsLepton(ParticleType particle) {
  if (IsChargedLepton(particle) || IsNeutrino(particle))
    return true;
  return false;
}
//_______________________________________________________________________
bool WimpSim::IsChargedLepton(ParticleType particle) {
  switch (particle) {
    case e_minus:
    case e_plus:
    case mu_minus:
    case mu_plus:
    case tau_minus:
    case tau_plus:
return true;
    default:
return false;
  }
}
//_______________________________________________________________________
bool WimpSim::IsNeutrino(ParticleType particle) {
  switch (particle) {
    case nu_e:
    case nu_e_bar:
    case nu_mu:
    case nu_mu_bar:
    case nu_tau:
    case nu_tau_bar:
return true;
    default:
return false;
  }
}
//_______________________________________________________________________
unsigned int WimpSim::SourceStringToInt(std::string sourcestr) {
  if (sourcestr == "Earth")
    return 2;
  else if (sourcestr == "Sun")
    return 1;
  else
    return 0;
}
  


//====================================================

I3WimpParams::~I3WimpParams() {}

I3WimpParams::I3WimpParams():
  source_(WimpSim::UNKNOWN),
  mass_(NAN),
  channel_(WimpSim::unknown),
  nu_weight_(NAN),
  lep_weight_(NAN),
  had_weight_(NAN),
  vgen_(NAN),
  time_(I3Time()) {}

I3WimpParams::I3WimpParams(WimpSim::SourceType source,
                           double mass,
                           WimpSim::DecayChannel channel,
                           double nu_weight,
                           double lep_weight,
                           double had_weight,
                           double vgen,
                           double aproj,
                           I3Time time):
  source_(source),
  mass_(mass),
  channel_(channel),
  nu_weight_(nu_weight),
  lep_weight_(lep_weight),
  had_weight_(had_weight),
  vgen_(vgen),
  aproj_(aproj),
  time_(time) {}

I3WimpParams::I3WimpParams(std::string source,
                           double mass,
                           int channel,
                           double nu_weight,
                           double lep_weight,
                           double had_weight,
                           double vgen,
                           double aproj,
                           I3Time time):
  mass_(mass),
  nu_weight_(nu_weight),
  lep_weight_(lep_weight),
  had_weight_(had_weight),
  vgen_(vgen),
  aproj_(aproj),
  time_(time)
{
  boost::to_upper(source);
  if (source=="SUN")
    source_ = WimpSim::SUN;
  else if (source=="EARTH")
    source_ = WimpSim::EARTH;
  else
    source_ = WimpSim::UNKNOWN;
  channel_=WimpSim::ConvertToDecayChannel(channel);
}

WimpSim::SourceType I3WimpParams::GetSource() const {return source_;}

std::string I3WimpParams::GetSourceString() const {
  if (source_ == WimpSim::EARTH)
    return "Earth";
  if (source_ == WimpSim::SUN)
    return "Sun";
  if (source_ == WimpSim::UNKNOWN)
    return "UNKNOWN";
  else
    return "";
}

double I3WimpParams::GetMass() const {return mass_;}

WimpSim::DecayChannel I3WimpParams::GetChannel() const {return channel_;}

int I3WimpParams::GetChannelInt () const {return (channel_);}
double I3WimpParams::GetNuWeight() const {return nu_weight_;}
double I3WimpParams::GetLepWeight() const {return lep_weight_;}
double I3WimpParams::GetHadWeight() const {return had_weight_;}

double I3WimpParams::GetLepWeightE47() const {return lep_weight_*1E47;}
double I3WimpParams::GetHadWeightE47() const {return had_weight_*1E47;}

double I3WimpParams::GetVgen() const {return vgen_;}

double I3WimpParams::GetAproj() const {return aproj_;}

I3Time I3WimpParams::GetTime() const {return time_;}

void I3WimpParams::SetSource(const WimpSim::SourceType source) {source_ = source;}
void I3WimpParams::SetSource(std::string source){
  boost::to_upper(source);
  if (source == "EARTH")
    source_= WimpSim::EARTH;
  else if (source == "SUN")
    source_= WimpSim::SUN;
  else if (source == "UNKNOWN")
    source_= WimpSim::UNKNOWN;
  else
    source_= WimpSim::UNKNOWN;
}
void I3WimpParams::SetSource(unsigned int source) {
  if (source == 2)
    source_= WimpSim::EARTH;
  else if (source == 1)
    source_= WimpSim::SUN;
  else if (source == 0)
    source_= WimpSim::UNKNOWN;
  else
    source_= WimpSim::UNKNOWN;
}


void I3WimpParams::SetMass(const double mass) {mass_ = mass;}

void I3WimpParams::SetChannel(const WimpSim::DecayChannel channel) {channel_ = channel;}// should go in unification with I3Particle

void I3WimpParams::SetNuWeight(const double nu_weight) {nu_weight_ = nu_weight;}
void I3WimpParams::SetLepWeight(const double lep_weight) {lep_weight_ = lep_weight;}
void I3WimpParams::SetHadWeight(const double had_weight) {had_weight_ = had_weight;}

void I3WimpParams::SetLepWeightE47(const double lep_weight) {lep_weight_ = lep_weight*1E-47;}
void I3WimpParams::SetHadWeightE47(const double had_weight) {had_weight_ = had_weight*1E-47;}

void I3WimpParams::SetVgen(const double vgen) {vgen_ = vgen;}

void I3WimpParams::SetAproj(const double aproj) {aproj_ = aproj;}

void I3WimpParams::SetTime(const I3Time time) {time_ = time;}
void I3WimpParams::SetTime(const double mjd) {time_ = I3Time(mjd);} //use I3Time.Constructor

template <class Archive>
void I3WimpParams::serialize(Archive& ar, unsigned version) {
  if (version>i3wimpparams_version_)
    log_fatal("Attempting to read version %u from file but running version %u of I3WimpParams class.",version,i3wimpparams_version_);

  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("Source", source_);
  ar & make_nvp("Mass",  mass_);
  ar & make_nvp("Channel", channel_);
  ar & make_nvp("Nu_Weight", nu_weight_);
  ar & make_nvp("Lep_Weight", lep_weight_);
  ar & make_nvp("Had_Weight", had_weight_);
  ar & make_nvp("Vgen", vgen_);
  ar & make_nvp("Aproj", aproj_);
  ar & make_nvp("Time", time_);
}

//-----------------------------------------------------------

std::ostream& operator<<(std::ostream& oss, const I3WimpParams& d){
  return(d.Print(oss));
}

std::ostream& I3WimpParams::Print(std::ostream& oss) const{
  oss << "[ I3WimpParams ::" << std::endl
      << "     Source : " << GetSourceString() << std::endl
      << " Mass (GeV) : " << GetMass()/I3Units::GeV << std::endl
      << "    Channel : " << GetChannel() << std::endl
      << "  Nu_Weight : " << GetNuWeight() << std::endl
      << " Lep_Weight : " << GetLepWeight() << std::endl
      << " Had_Weight : " << GetHadWeight() << std::endl
      << "       Vgen : " << GetVgen() << std::endl
      << "      Aproj : " << GetAproj() << std::endl
      << "       Time : " << GetTime().GetUTCString() << std::endl
      << "]";
  return oss;
}

I3_SERIALIZABLE(I3WimpParams);
