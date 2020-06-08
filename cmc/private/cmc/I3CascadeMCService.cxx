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
#include "cmc/I3CascadeSplit.h"
#include "cmc/I3CascadeMuonSplit.h"
#include <cmc/I3CascadeMCService.h>
#include "I3CascadeMCCommon.h"
// icetray includes
#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include <boost/foreach.hpp>

using I3CascadeMCCommon::IsCascade;
using I3CascadeMCCommon::IsHadron;
using I3CascadeMCCommon::HadronEnergyCorrection;

using namespace std;

/**
 * Constructor
 */
I3CascadeMCService::I3CascadeMCService(I3RandomServicePtr r) : 
  random_(r),
  splitter_( new I3CascadeSplit( r ) ),
  muonSplitter_( new I3CascadeMuonSplit( r ) )
 {
  log_trace("Entering Constructor");
  // create splitter
  SetStepWidth(I3CascadeSplit::DEFAULT_STEP_WIDTH);
  SetEnergyThresholdSimulation(1 * I3Units::PeV);
  splitter_->SetSegmentMaxEnergy(1 * I3Units::PeV);

  // create muon splitter
  SetThresholdSplit(I3CascadeParametrization::DEFAULT_THRESHOLD);
  SetEnergyThresholdMuons(I3CascadeMuonSplit::DEFAULT_ENERGY_THRESHOLD);
  SetMaxMuons(I3CascadeMuonSplit::DEFAULT_MAX_MUONS);
}

// interface defined by I3PropagatorService
std::vector<I3Particle> I3CascadeMCService::Propagate(I3Particle& p, DiagnosticMapPtr frame, I3FramePtr) {
  std::vector<I3Particle> retval;
  Simulate(p, retval);
  return retval;
}

// interface defined by I3PropagatorService
void I3CascadeMCService::SetRandomNumberGenerator(I3RandomServicePtr random) {
  random_ = random;
  splitter_->SetRandomNumberGenerator(random);
  muonSplitter_->SetRandomNumberGenerator(random);
}

void I3CascadeMCService::SetEnergyThresholdMuons( double t){ 
  energyThresholdMuons_ = t;
  muonSplitter_->SetEnergyThresholdMuons(energyThresholdMuons_);
};

void I3CascadeMCService::SetMaxMuons( int i ){ 
  maxMuons_ = i;
  muonSplitter_->SetMaxMuons(maxMuons_);
};

void I3CascadeMCService::SetThresholdSplit(double t ){ 
  energyThresholdSplit_ = t; 
};

void I3CascadeMCService::SetEnergyThresholdSimulation(double t){ 
  energyThresholdSimulation_ = t;
  splitter_->SetSimulationThreshold(energyThresholdSimulation_);
};

void I3CascadeMCService::SetStepWidth(int w ){ 
  stepWidth_ = w;
  splitter_->SetStepWidth(stepWidth_);
};

// replace the cascade in the given tree pointed to by iter with a list of sub-cascades
// returns the iterator pointing to the last added particle
void I3CascadeMCService::Simulate(I3Particle& cascade, vector<I3Particle>& daughters) {

  
  if (!IsCascade(cascade))
    log_fatal("%s is not a cascade!", cascade.GetTypeString().c_str());

  cascade.SetLength(0.);
  I3Particle new_cascade = cascade.Clone();

  // if particle is a hadron shower rescale hadron energy and 
  // substitue hadron shower with EM shower
  if ( IsHadron(new_cascade) ) {
    log_debug("Got a %.1e GeV hadronic shower (%s)", new_cascade.GetEnergy(), new_cascade.GetTypeString().c_str());
    
    if (maxMuons_ > 0 && 
	new_cascade.GetEnergy() >= energyThresholdMuons_ ) {
      // and the cascade energy is sufficient
      // and requested number of muons is > 0
      vector<I3Particle> cascadeMuons 
	= muonSplitter_->GenerateMuons(new_cascade);
      log_trace("Generate %zi muons", cascadeMuons.size());
      BOOST_FOREACH(I3Particle& c, cascadeMuons) daughters.push_back(c);
    }
    
    // don't simulate hadronic shower profile, it's not affected by the LPM effect
    splitter_->EnableSimulation(false);

  } else {
    // always simulate em-shower (if the energy is above the configured threshold)
    splitter_->EnableSimulation(true);
  }

  // split the cascade into a list of sub-cascades
  if (cascade.GetEnergy() > energyThresholdSplit_) {
    cascade.SetShape(I3Particle::Dark);
    vector<I3Particle> subCascades = splitter_->SplitCascade(new_cascade);
    log_debug("Adding %zi sub-cascades", subCascades.size());
    BOOST_FOREACH(I3Particle& c, subCascades) daughters.push_back(c);
  // below the threshold, only emit a daughter if the visible energy had to be adjusted
  } else if ((cascade.GetType() != I3Particle::EMinus) &&
    (new_cascade.GetEnergy() != cascade.GetEnergy())) {
    cascade.SetShape(I3Particle::Dark);
    daughters.push_back(new_cascade);
  }
  
  log_debug("Added a total of %zu daugthers", daughters.size());

}


