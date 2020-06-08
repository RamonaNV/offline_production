/**
 * Cascade longitudinal energy profile parametrization
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeParametrization.h 32874 2007-06-04 18:22:36Z bvoigt $
 *
 * @version $Revision: 32874 $
 * @date $LastChangedDate: 2007-06-04 12:22:36 -0600 (Mon, 04 Jun 2007) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: bvoigt $
 */

#ifndef I3_CASCADE_PARAMETRIZATION_H_INCLUDED
#define I3_CASCADE_PARAMETRIZATION_H_INCLUDED

//Standard C/C++ includes
#include <vector>

// local includes
#include "cmc/I3CascadeDevelopment.h"
#include "cmc/I3CascadeMCCommon.h"

/**
 *@brief This class implements the longitudinal shower profile of electromagnetic cascades
 *       for energies below the LPM Threshold (roughly 1PeV in Ice)
 *
 * The parametrization is take from the Particle Data Booklet, Section 26.5
 * The parameters a,b of the dE/dt function described their have been determined by
 * Christoher Wiebusch using GEANT 3. They are taken from his thesis.
 *
 *
 * @author Bernhard Voigt
 */
class I3CascadeParametrization : public I3CascadeDevelopment, public I3CascadeMCCommon::Sampler {

  // define the logger name
  SET_LOGGER("I3CascadeParametrization");

 public:

  /**
   * Default constructor
   */
  I3CascadeParametrization();

  /**
   * Default destructor
   */
  ~I3CascadeParametrization() {}

  /**
   * Evaluates the longitudinal shower profile parametrization for an EM shower of given energy
   *
   * The parametrization is taken from PDG Booklet Section (Particle Passage through matter)
   *
   * \f$\frac{dE}{dt} = E\ b\ \frac{(bt)^{a-1}\ e^{-bt}}{\Gamma(a)}\f$
   *
   * The parameters used are a, b which have been determined by Christopher Wiebusch and quoted
   * in his thesis.
   *
   * The energy deposit is calculated in steps of one radiation length. It is internaly stored
   * and can be accessed using the I3CascadeDevelopment::GetEnergyDeposit method.
   *
   * @todo add shower type parameter (EM or Hardonic)
   *
   * @param energy energy of the shower
   */
  void Simulate(I3Particle::ParticleType type, double energy);

  /**
   * Calculates the energy loss function dE/dt
   *
   * @param energy of the shower
   * @param t is the depth in radiation length units
   */
  double dE_dt(double energy, double t);

  /**
   * Sets the threshold for the parametrization
   *
   * Calculation of energy deposit is quit, when the
   * deposit in one radiation length step is less than 
   * this threshold, default is 1 TeV
   *
   * @param threshold
   */
  void SetThreshold(double threshold) {
    threshold_ = threshold;
  }

  /**
   * Returns the threshold of energy deposit calculation
   */
  double GetThreshold() {
    return threshold_;
  }

  /**
   * @brief threshold for energy deposit calculation
   */
  static const double DEFAULT_THRESHOLD;


 private:

  /**
   * @brief Threshold for calculation energy deposit
   *
   * Calculation of energy deposit is quit, when the
   * deposit in one radiation length step is less than 
   * this threshold
   */
  double threshold_;

  /**
   * @brief parameter for 'a' parameter calculation, no units.
   *
   * From Wiebusch thesis
   */
  static const double A_0;
  static const double A_1;

  /**
   * @brief paramter for dE/dt function
   *
   * From Wiebusch thesis
   */
  static const double B;
};
#endif
