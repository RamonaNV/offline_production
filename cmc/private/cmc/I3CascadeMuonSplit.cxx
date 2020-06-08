/**
 * brief Implementation I3CascadeMuonSplit class
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSplit.cxx 33092 2007-06-11 14:49:59Z bvoigt $
 *
 * @version $Revision: 33092 $
 * @date $LastChangedDate: 2007-06-11 16:49:59 +0200 (Mon, 11 Jun 2007) $
 * @author Sebastian Panknin <panknin@physik.hu-berlin.de>  
 * Last changed by: $LastChangedBy: bvoigt $
 */

//Standard C/C++ includes
#include <cmath>
#include <vector>
#include <algorithm>

// local includes
#include "I3CascadeMuonSplit.h"

//IceTray includes
#include "dataclasses/I3Constants.h"
#include "icetray/I3Units.h"
#include "phys-services/I3GSLRandomService.h"

using namespace std;

// Class constants
const double I3CascadeMuonSplit::DEFAULT_ENERGY_THRESHOLD 
  = 1 * I3Units::GeV;
const int I3CascadeMuonSplit::DEFAULT_MAX_MUONS 
  = 10;  

I3CascadeMuonSplit::I3CascadeMuonSplit(I3RandomServicePtr random):
  I3CascadeMCCommon::Sampler(random),
  energyThresholdMuons_(DEFAULT_ENERGY_THRESHOLD),
  maxMuons_(DEFAULT_MAX_MUONS)
{
}

I3CascadeMuonSplit::~I3CascadeMuonSplit() {}

void I3CascadeMuonSplit::SetEnergyThresholdMuons(double threshold) {
  energyThresholdMuons_ = threshold;
}

void I3CascadeMuonSplit::SetMaxMuons(int max) {
  maxMuons_ = max;
}


// Implementation of the cascade muon split method.
vector<I3Particle> I3CascadeMuonSplit::GenerateMuons(I3Particle &cascade) {
  log_trace("Entering GenerateMuons");

  vector<I3Particle> particleList;

  if( !std::isfinite(cascade.GetEnergy()) ){
    log_error("The cascade has a pathological energy = %e",
	      cascade.GetEnergy());
    log_error("Skipping this cascade.");
    return particleList;
  }

  // Parameters for the parametresation according to Corsica data
  // from Sebastian Panknin and Marek Kowalski, 2007
  double nint_muon = 0.358/I3Units::TeV; 
  double gammint_muon = 1.74; 

  // calculate the mean number of muons expected in this cascade
  double meanMuonNumber =  cascade.GetEnergy() * nint_muon 
    * pow(energyThresholdMuons_/I3Units::GeV, -gammint_muon);
  log_trace("Mean number of muons: %f", meanMuonNumber);

  // muon production follows poisson statistic
  int numberOfMuons = random_->Poisson(meanMuonNumber);
  log_debug("Number of muons: %i", numberOfMuons);


  // Calculate an energy for each muon
  vector<double> muonEnergies;
  double totalEnergy = 0.;
  for (int i = 0; i < numberOfMuons; i++) {
    double muonEnergy = pow(10.,
                            log10(energyThresholdMuons_) 
			    - log(random_->Uniform())
			      /(gammint_muon*log(10)));
    log_trace("%i Muon Energy: %.1f", i, muonEnergy); 

    // if the sum of all energies exceeds the cascade energy
    // assign only the energy which is still available
    if (totalEnergy + muonEnergy > cascade.GetEnergy()) {
      muonEnergy = cascade.GetEnergy() - totalEnergy;
    }
    muonEnergies.push_back(muonEnergy);
    totalEnergy += muonEnergy;
  }

  // sort the energies in reverse order
  sort(muonEnergies.begin(),muonEnergies.end()); // std::algorithm
  reverse(muonEnergies.begin(), muonEnergies.end()); // std::algorithm

  // only generate the maximum requested number of muons
  if (numberOfMuons > maxMuons_) {
    log_trace("Number of muons exceeds maxMuons %i option, "\
	      "maxMuons will be used", maxMuons_);
    numberOfMuons = maxMuons_;
  }

  // reset the total energy 
  // - only numberOfMuons should be taken into account
  totalEnergy = 0;

  for (int i=0; i < numberOfMuons; i++) {

    // copy the original to keep the direction and origin
    I3Particle muon; 
    // set the muon properties
    muon.SetLocationType(cascade.GetLocationType()); 
    muon.SetDir(cascade.GetDir());
    muon.SetPos(cascade.GetPos());
    muon.SetTime(cascade.GetTime());
    muon.SetType(I3Particle::MuPlus);
    muon.SetEnergy(muonEnergies[i]);
    //the muons now need to be propagated mmc
    muon.SetLength(NAN);
    particleList.push_back(muon);
    log_debug("Added Muon - E: %.1f Length: %f", 
	      muon.GetEnergy(), muon.GetLength());

    // accumulate the energy that went into muons
    totalEnergy += muonEnergies[i];
  }
  log_debug("Added %i muons with a total energy of %.1f to the event", 
            numberOfMuons, totalEnergy);

  // substract muon energy from the cascade
  cascade.SetEnergy(cascade.GetEnergy() - totalEnergy);
  log_debug("Muon Energy substracted - Particle E: %.1f", 
	    cascade.GetEnergy());

  log_trace("Leaving GenerateMuons");
  return particleList;
}




    
    

  


