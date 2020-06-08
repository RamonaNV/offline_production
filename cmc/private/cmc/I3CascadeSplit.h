/**
 * Cascade splitting for a realistic treatment of the longitudinal development of cascasdes
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSplit.h 77043 2011-06-21 19:43:30Z olivas $
 *
 * @version $Revision: 77043 $
 * @date $LastChangedDate: 2011-06-21 13:43:30 -0600 (Tue, 21 Jun 2011) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: olivas $
 */

#ifndef I3CASCADESPLIT_h_
#define I3CASCADESPLIT_h_

//Standard C/C++ includes
#include <cmath>
#include <vector>

// Local includes
#include "cmc/I3CascadeDevelopment.h"
#include "cmc/I3CascadeParametrization.h"
#include "cmc/I3CascadeSimulation.h"
#include "cmc/I3CascadeMCCommon.h"

// icetray includes
#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"

/**
 *@brief This class defines methods to split a cascade into sub-showers
 *
 * There is only one method of interest to the user which is the split method.
 * It takes a I3Particle of a cascade typ, splits it into 
 * a set of sub-cascades and returns the corresponding I3Particles.
 *
 * The parametrizations used for the splitting are defined in a I3CascadeParametrization object
 * which is passed to the constructor of this class. Update docu!!!
 *
 * @author Bernhard Voigt
 */
class I3CascadeSplit : public I3CascadeMCCommon::Sampler {

  // set logging name for icetray
  SET_LOGGER("I3CascadeSplit");

 public:

  /**
   * Creates a I3CascadeSplit object
   *
   * @param random service pointer
   */
  I3CascadeSplit(I3RandomServicePtr random);

  /**
   * Copy constructor 
   *
   * @param cascadeSplit I3CascadeSplit
   */
private:
  I3CascadeSplit(const I3CascadeSplit& cascadeSplit);
public:
  virtual void SetRandomNumberGenerator(I3RandomServicePtr r);
  
  /**
   * Default destructor
   * 
   */
  virtual ~I3CascadeSplit();
  
  /**
   * @brief Splits a given I3Particle into multiple I3Partitcle.
   *
   * A cascade-like particle is split into multiple sub-cascades to 
   * simulate the longitudinal development of the shower.
   * The length and shower profile are obtained from the I3CascadeDevelopment object
   * passed to the constructor of this object.
   * A vector of I3Particles is returned, with individual particle information, each of
   * a EM cascade type (I3Particle::EMinus). 
   *
   * @param cascade I3Particle
   * @return vector of I3Particles 
   */
  std::vector<I3Particle> SplitCascade(I3Particle &cascade);

  /**
   * Set the step width between sub-cascades
   *
   * @param stepWidth distance between sub-cascades in units of radiation length
   */
  void SetStepWidth(int stepWidth) {
    stepWidth_ = stepWidth;
  }

  /**
   * Returns the step width, distance between sub-cascades
   */
  int GetStepWidth() {
    return stepWidth_;
  }

  /**
   * Simulate high energy cascades, rather than using
   * the enrgy loss profile parameterisation
   *
   * @param state
   */
  void EnableSimulation(bool state) {
    enableSimulation_ = state;
  }

  /**
   * Returns true when the simulation for high energy particles
   * is enabled
   */
  bool IsEnableSimulation() {
    return enableSimulation_;
  }

  /**
   * Simulate cascades with energies above this threshold
   *
   * @param threshold energy threshold above which cascades are simulated
   */
  void SetSimulationThreshold(double threshold) {
    simulationThreshold_ = threshold;
  }

  /**
   * Returns the threshold above which cascades are simulated
   */
  double GetSimulationThreshold() {
    return simulationThreshold_;
  }

  /**
   * Ensure that all output cascades are split up to be no larger than this
   *
   * @param threshold energy threshold above which ouptut cascades are split
   */
  void SetSegmentMaxEnergy(double threshold) {
    segmentMaxEnergy_ = threshold;
  }
  
  /**
   * Returns the threshold above which ouptut cascades are split
   */
  double GetSegmentMaxEnergy() {
    return segmentMaxEnergy_;
  }

  static const double DEFAULT_SIMULATION_THRESHOLD;
  static const int DEFAULT_STEP_WIDTH;
  static const double DEFAULT_MAX_ENERGY;

 private:

  /**
   * Creates a I3Particle from a given Particle with given energy, length in a distance of 
   * steps * radiation length away from the original I3Particle position in the direction 
   * of the original I3Particle
   *
   * @param particle 
   * @param energy
   * @param steps
   * @param length
   * @return newly created I3Particle
   */
  I3Particle CreateParticle(I3Particle &particle, double energy, int steps, int length);

  /**
   * @brief PSI_CascadeDevelopmentParametrization object providing a parametrization of 
   * the longitudinal shower profile and shower length
   */
  I3CascadeParametrization cascadeParametrization_;

  /**
   * @brief PSI_CascadeDevelopment object to obtain the longitudinal shower profile and 
   * length from a simulation
   */
  I3CascadeSimulation cascadeSimulation_;

  /**
   * @brief switch to enable the cascade simulation
   */
  bool enableSimulation_;

  /**
   * @brief threshold energy for cascade simulation, if simulation is enable showers above 
   * this energy will be simulated
   */
  double simulationThreshold_;

  /**
   * @brief distance between sub-cascades in units of radiation length
   */
  int stepWidth_;
  
  /**
   * @brief maximum energy to output in a single cascade
   *
   * Output cascades which would exceed this energy will be further split into
   * subcascades of this energy, all at the same position.
   */
  double segmentMaxEnergy_;

};
#endif
