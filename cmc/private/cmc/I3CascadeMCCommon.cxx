/**
 *@brief Implementation I3CascadeMCService class
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeMCService.cxx 66313 2010-08-23 16:09:43Z olivas $
 *
 * @version $Revision: 66313 $
 * @date $LastChangedDate: 2010-08-23 12:09:43 -0400 (Mon, 23 Aug 2010) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: olivas $
 */

// local includes
#include "I3CascadeMCCommon.h"

// icetray includes
#include "icetray/I3Units.h"
#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"

// check whether the particle is of any type which produces a shower
bool I3CascadeMCCommon::IsCascade(const I3Particle& particle) {
  log_trace("Particle type: %s shape: %s Energy: %.3e", 
            particle.GetTypeString().c_str(), 
            particle.GetShapeString().c_str(),
            particle.GetEnergy());

  // check whether the particle is of any type which produces a shower
  // I don't use I3Particle::IsCascade since it doesn't include all the types
  // checked here
  if ( particle.GetType() == I3Particle::Gamma ||
       particle.GetType() == I3Particle::EPlus ||
       particle.GetType() == I3Particle::EMinus ||
       particle.GetType() == I3Particle::Brems ||
       particle.GetType() == I3Particle::DeltaE ||
       particle.GetType() == I3Particle::PairProd ||
       particle.GetType() == I3Particle::NuclInt || 
       particle.GetType() == I3Particle::Hadrons ||
       particle.GetType() == I3Particle::PPlus ||
       particle.GetType() == I3Particle::PMinus ||
       particle.GetType() == I3Particle::Pi0 ||
       particle.GetType() == I3Particle::PiPlus ||
       particle.GetType() == I3Particle::PiMinus ) {

    log_trace("Particle is a casacade!");
    return true;
  } else {
    log_trace("Particle is not a casacade!");
    return false;
  }
 }

// checks whether the particle is a hadronic shower
bool I3CascadeMCCommon::IsHadron(const I3Particle& particle) {
  if( particle.GetType() == I3Particle::Hadrons ||
      particle.GetType()==I3Particle::PiPlus ||
      particle.GetType()==I3Particle::PiMinus ||
      particle.GetType()==I3Particle::NuclInt ||
      particle.GetType() == I3Particle::PPlus ||
      particle.GetType() == I3Particle::PMinus ) {

    return true;
  } else {
    return false;
  }
}

//Scales energy of hadronic shower to fraction of EM shower according to formula by M.Kowalski
double I3CascadeMCCommon::HadronEnergyCorrection(double energy, I3RandomServicePtr rand) {
  // numbers are taken from M.Kowalski's internal report
  
  // hadron energy scaling factor
  // this is how much less light is produced by a hadron shower, compared to an EM shower
  // Based on looking at figure 4.4 in M. Kowalski's AMANDA Internal Report 
  // (found here http://icecube.berkeley.edu/amanda-private/reports/20020803-track.pdf)
  // Although the equations in the report are for ln (log), calculation reveals the 
  // points in figure 4.4 (and 3.2) are best fit by log10
  // Also, equation 3.4 should be RMS=F*RMS0*(log10(E))**-gamma according to this e-mail 
  // from Marek: http://lists.icecube.wisc.edu/pipermail/ice3cascade/2010-July/002588.html
  energy = std::max(10.0, energy);
  double energyScalingFactor = 1 + pow( (energy/0.399), -0.130 )*( 0.467 - 1 );
  
  if( !std::isnormal(energy) ){
    log_error("energy (%e) is not normal...returning the same energy", energy);
    return energy;
  }

  double logenergy = log10(energy);
  // width of the scaling factor distribution
  double std = energyScalingFactor * 0.379 * pow(logenergy, -1.160);
  // smear scaling factor
  // we don't want numbers that are below zero or greater than 1, so keep sampling if it's out of range
  double newenergyScalingFactor=-1;
  while( newenergyScalingFactor < 0
	 || newenergyScalingFactor > 1 ){
    newenergyScalingFactor = rand->Gaus(energyScalingFactor, std);
  }
  energyScalingFactor=newenergyScalingFactor;

  // scale energy
  energy *= energyScalingFactor;

  return energy;
}


