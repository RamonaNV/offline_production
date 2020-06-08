/**
 * Cascade longitudinal base class for parametrization and simulation classes
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeDevelopment.h 32874 2007-06-04 18:22:36Z bvoigt $
 *
 * @version $Revision: 32874 $
 * @date $LastChangedDate: 2007-06-04 12:22:36 -0600 (Mon, 04 Jun 2007) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: bvoigt $
 */

#ifndef I3_CASCADE_DEVELOPMENT_H_INCLUDED
#define I3_CASCADE_DEVELOPMENT_H_INCLUDED

// icetray includes
#include "icetray/I3Logging.h"
#include "dataclasses/physics/I3Particle.h"

/**
 *@brief This abstract class defines the interface of a cascade development description
 *
 *
 * An implementation of this class would define a Simulate method that calculates the
 * length and the longitudinal energy deposit profile of a cascade of given energy.
 *
 * After this method has been invoked, the length and longitudinal profile are accessible
 * via the corresponding methods of this class.
 *
 * Implementations can be either a pure parametrization or a simulation of the shower development.
 *
 * @author Bernhard Voigt
 */
class I3CascadeDevelopment {

  SET_LOGGER("I3CascadeDevelopment");

 public:

  /**
   * Virtual destructor
   */
  virtual ~I3CascadeDevelopment() {}

  /**
   * Simulates the cascade development
   *
   * This method must be called in order to get the length and/or
   * longitudinal profile.
   *
   * @todo Add shower type parameter (EM or Hardonic)
   *
   * @param energy energy of the shower
   */
  virtual void Simulate(I3Particle::ParticleType type, double energy) = 0;

  /**
   * Returns the lenght of a shower units of radiation length
   * 
   * I3CascadeDevelopment::Simulate must be called previously in order to get the shower length
   *
   * @return the length of the shower in units of radiation length
   */
  double GetLength();

    /**
     * Returns the energy loss of the shower between a and b (units are radiation length)
     * 
     * I3CascadeDevelopment::Simualte must be called beforehand.
     *
     * @param a start position (units are radiation length)
     * @param b end position (units are radiation length)
     * @return fractional energy loss between a and b (excluding b)
     */
  double GetEnergyDeposit(int a, int b);

  /**
   * @brief Returns the energy loss profile filled in the simulate method
   *
   */
  std::vector<double> GetEnergyLossProfile() {
    return energyLossProfile_;
  }

  /**
   * @brief radiation length in ice in meter
   */
  static const double RADIATION_LENGTH; 

 protected:

  /**
   * @brief pointer to the vector holding the energy loss profile at every radiation length
   *        filled in Simulate method.
   */
  std::vector<double> energyLossProfile_;

};
#endif
