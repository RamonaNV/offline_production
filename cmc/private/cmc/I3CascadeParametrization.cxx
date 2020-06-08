/**
 *@brief Implementation I3CascadeParametrization class
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeParametrization.cxx 150701 2016-10-12 14:22:25Z jvansanten $
 *
 * @version $Revision: 150701 $
 * @date $LastChangedDate: 2016-10-12 08:22:25 -0600 (Wed, 12 Oct 2016) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: jvansanten $
 */

//Standard C/C++ includes
#include <cmath>
#include <vector>
#include <numeric>

// GNU Scientific Library includes
#include <gsl/gsl_sf_gamma.h>  // gamma function from special-function package

// icetray includes 
#include "icetray/I3Units.h"
#include "sim-services/I3SimConstants.h"
#include "phys-services/I3RandomService.h"

// local includes
#include "I3CascadeDevelopment.h"
#include "I3CascadeParametrization.h"

// definition of class constants
const double I3CascadeParametrization::DEFAULT_THRESHOLD = 1 * I3Units::TeV;

// Default constructor
I3CascadeParametrization::I3CascadeParametrization() :
  I3CascadeDevelopment(),
  threshold_(I3CascadeParametrization::DEFAULT_THRESHOLD) {
}

// Implementation of the simulate method that evaluates the energy loss parametrization
void I3CascadeParametrization::Simulate(I3Particle::ParticleType type, double energy) {
  /* compute energy loss function for every radiation length
   * the parametrization extends to infinity, leaving very small energy fractions
   * in the tail of the cascade
   * thus compute the energy loss to a certain point and distribute the missing
   * energy along the shower
   */
  
  // init a fresh energy profile
  log_trace("Deleting old profile");
  energyLossProfile_.clear();
  double energyLoss = 0.0;
  double lastEnergyLoss = 0.0;
  double remainingEnergy = energy;
  double step = 1;  // start at a depth of one radiation length

  // energy loss formula taken from PDB booklet (Sec. 26.5)
  // dE/dt = E * b * ( (bt)^(a-1) * exp(-bt) )/Gamma(a)
  I3SimConstants::ShowerParameters params(type, energy);
  double a = params.a;
  double b = I3CascadeDevelopment::RADIATION_LENGTH/params.b;
  
  double f = 1.;
  if (params.emScaleSigma != 0.) {
      do {
          f=params.emScale +
              params.emScaleSigma*random_->Gaus(0.,1.);
      } while((f<0.) || (1.<f));
  }

  if (!(a > 0 && a < GSL_SF_GAMMA_XMAX)) {
    log_error("Gamma function argument a=%f out of range!", a);
    return;
  }
  const double norm = energy*b/gsl_sf_gamma(a);
  
  // add energyLoss to the energyLossProfile_ vector
  // until energy loss is less than threshold except, this is still
  // the raising edge of the parametrization
  do {
    // keep last energy loss to see whether this is still the raising edge
    lastEnergyLoss = energyLoss;
    energyLoss = norm * pow(b*step, a-1) * exp(-b*step);
    energyLossProfile_.push_back(f*energyLoss);
    // this energyLoss has been taken into account
    remainingEnergy -= energyLoss;
    step++; // next step (unit is radiation length)
    log_trace("Missing energy %.3e, threshold is %.3e", remainingEnergy, threshold_);
  } while (lastEnergyLoss < energyLoss || energyLoss > threshold_);

  // distribute the missing energy along the shower
  std::vector<double>::iterator iter = energyLossProfile_.begin();
  double sum = accumulate(energyLossProfile_.begin(), energyLossProfile_.end(), 0.0);
  for (; iter != energyLossProfile_.end(); iter++) {
    // add fraction of remaining energy 
    *iter += (*iter/sum) * f*remainingEnergy;
  }
  log_debug("Length: %zi Sum: %.3e", energyLossProfile_.size(), 
	    accumulate(energyLossProfile_.begin(), 
		       energyLossProfile_.end(), 0.0));
}

// implementaion of the energy loss profile
// formula taken from PDB booklet (Sec. 26.5)
// parameters defined in header file (taken from Wiebusch Thesis)
// dE/dt = E * b * ( (bt)^(a-1) * exp(-bt) )/Gamma(a)
double I3CascadeParametrization::dE_dt(double energy, double t) {
  I3SimConstants::ShowerParameters params(I3Particle::EMinus, energy);
  double a = params.a;
  double b = I3CascadeDevelopment::RADIATION_LENGTH/params.b;
  log_debug("Energy %.3e - parameter to gamma function %f", energy, a);
  if (!(a > 0 && a < GSL_SF_GAMMA_XMAX)) {
    log_error("Parameter to gamma function out of range a=%f."
	     " Energy might be too low or too large e=%e", a, energy);
    return 0.0;
  }
  double norm = 1/gsl_sf_gamma(a);
  return norm * energy * b * pow(b*t, a-1) * exp(-b*t); 
}



  

