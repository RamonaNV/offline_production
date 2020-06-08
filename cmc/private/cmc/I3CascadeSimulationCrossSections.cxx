/**
 *@brief Implementation I3CascadeSimulationCrossSections static methods
 *
 * All formulars are taken from S.Klein's review. There's not much documentation
 * since it's just copied from the paper and there's no logic, just formulars.
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSimulationCrossSections.cxx 148463 2016-07-15 06:30:57Z gary.binder $
 *
 * @version $Revision: 148463 $
 * @date $LastChangedDate: 2016-07-15 00:30:57 -0600 (Fri, 15 Jul 2016) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: gary.binder $
 */

#include <cmath>

// local includes
#include "I3CascadeSimulationCrossSections.h"

// icetray includes
#include "dataclasses/I3Constants.h"

// class constants
const int I3CascadeSimulationCrossSections::NUMBER_OF_COMPONENTS = 3;
const int I3CascadeSimulationCrossSections::Z[I3CascadeSimulationCrossSections::NUMBER_OF_COMPONENTS] = {1, 1, 8};
const double I3CascadeSimulationCrossSections::X0 = 36.0; // Unit is cm
const double I3CascadeSimulationCrossSections::ALPHA = 1/137.;
const double I3CascadeSimulationCrossSections::ELECTRON_RADIUS = 2.817e-13; // Unit is cm
const double I3CascadeSimulationCrossSections::N = 6.022e23; // Unit is 1/mol
const double I3CascadeSimulationCrossSections::DENSITY = 0.917; // Unit g/cm**3
const double I3CascadeSimulationCrossSections::A = 18.0153; // Unit g/mol
const double I3CascadeSimulationCrossSections::Z_A_RATIO = .55509;

// number of nukleons, formular from Spencer's review
const double I3CascadeSimulationCrossSections::NUMBER_OF_NUKLEONS = 
  1/Z_A_RATIO * N * DENSITY/A * 4 * ALPHA * pow(ELECTRON_RADIUS, 2.0) * \
  MeanZ(2.0) * log(184.0/(MeanZ(1.0/3.0)));

const double I3CascadeSimulationCrossSections::ENERGY_THRESHOLD =  
  7.7e3 * 36.0 * I3Units::GeV; // this is 7.7e3 GeV * X0 [g/cm^2]


// calculate mean z: sum(components**power)/len(components)
double I3CascadeSimulationCrossSections::MeanZ(const double& exp){
  double mean = 0.0;
  for (int i=0; i < 3; i++) {
    mean = mean + pow(Z[i], exp);
  }
  //log_trace("Mean Z**%.2f: %.2f\n", exp, mean/NUMBER_OF_COMPONENTS);
  return mean/NUMBER_OF_COMPONENTS;
}  

double I3CascadeSimulationCrossSections::Xi(const double& s) {
  double s1 = MeanZ(2.0/3.0) * pow(184.0, 2.0);

  if (s < sqrt(2.0) * s1) {
    return 2.0;
  }

  if (sqrt(2.0) < s and s < 1) {
    double h = log(s) / log(sqrt(2.0)*s1);
    return 1 + h - (.08 * (1-h) * (1-pow(1-h,2.0))) / log(sqrt(2.0)*s1);
  }

  if (s >= 1 ) {
    return 1.0;
  }
  return .0;
}

double I3CascadeSimulationCrossSections::Phi(const double& s) {
  if (s <= .01) {
    return 6 * s * (1 - I3Constants::pi * s);
  }
   
  if (.01 < s and s < 2) {
    return 1 - exp(-6 * s * (1 + (3-I3Constants::pi) * s) + pow(s, 3.0) / \
                   (.623 + .796 * s + .658 * pow(s,2.0)));
  }

  if (s >= 2) {
    return 1.0;
  }
  return .0;
}

double I3CascadeSimulationCrossSections::Psi(const double& s) {
  if (s <= .01) {
    return 4 * s;
  }

  if (.01 < s and s < 2) {
    return 1 - exp(-4 * s - 8 * pow(s,2.0) / \
                   (1 + 3.96 * s + 4.97 * pow(s,2.0) - .05 * pow(s,3.0) + 7.5 * pow(s,4.0)));
  }

  if (s >= 1) {
    return 1.0;
  }
  return .0;
}

double I3CascadeSimulationCrossSections::G(const double& s) {
  return 3 * Psi(s) - 2 * Phi(s);
}

double I3CascadeSimulationCrossSections::BremsstrahlungCrossSection(const double& energy, const double& y) {
  // supression factor
  double s = sqrt(1.0/8.0 * ENERGY_THRESHOLD/energy * y/(1-y));

  // photon energy is the fraction of the primary energy
  double k = energy * y;

  return NUMBER_OF_NUKLEONS * energy * Xi(s) / (3.0*k) * \
    (pow(y,2.0) * G(s) + 2 * ((1 + pow(1-y,2.0)) * Phi(s)));
}

double I3CascadeSimulationCrossSections::PairProductionCrossSection(const double& energy, const double& y) {
  // supression factor
  double s = sqrt(1.0/8.0 * ENERGY_THRESHOLD/energy * 1/(y * (1-y)));

  return NUMBER_OF_NUKLEONS * Xi(s) / 3 * \
    (G(s) + 2 * (pow(y,2.0) + pow(1-y,2.0)) * Phi(s));
}
    

