/**
 * Cascade simulation based on pair production and bremsstrahlung interactions
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSimulation.h 35916 2007-08-23 09:34:12Z bvoigt $
 *
 * @version $Revision: 35916 $
 * @date $LastChangedDate: 2007-08-23 03:34:12 -0600 (Thu, 23 Aug 2007) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: bvoigt $
 */

#ifndef I3_CASCADE_SIMULATION_H_INCLUDED
#define I3_CASCADE_SIMULATION_H_INCLUDED

//Standard C/C++ includes
#include <vector>
#include <stack>

// Gnu Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

// local includes
#include "cmc/I3CascadeMCCommon.h"
#include "cmc/I3CascadeDevelopment.h"
#include "cmc/I3CascadeParametrization.h"
#include "cmc/I3MetropolisHastings.h"

// icetray includes
#include "phys-services/I3RandomService.h"

/**
 * @brief Local definition of a Particle struct
 *
 * The particle struct holds energy and postion informations
 * of particles in the simulation.
 */
struct Particle {
  static const short ELECTRON = 1;
  static const short PHOTON = 0;
  double energy;
  double x;
  short type;
  /**
   * Particle with energy, at position x (1-dimensional), of 
   * type type (1 is electron, 0 is photon)
   *
   * @param energy
   * @param x
   * @param type
   */
  Particle(double energy, double x, short type) : energy(energy), x(x), type(type) { }
};

/**
 * @brief Local definition of beta distribution
 *
 * The struct defines a pdf function, ie. the beta distribtuion
 * not normed, though.
 * It provides a rvs method to draw random samples from
 * the beta distribution using gsl_rng
 * 
 * The parameters used for the beta distribibution (alpha and beta)
 * are defined.
 *
 * Functions are passed to I3MetropolisHastings and need to be static.
 */
struct beta_dist {
  // define the logger name, use the same as for the simulation class
  SET_LOGGER("I3CascadeSimulation");

  // parameter to beta function
  static double alpha;
  static double beta;

  /**
   * Beta Function \f$ x^{\alpha-1} (1-x)^{\beta-1}\f$
   *
   * @param x
   */
  static double pdf(const double& x) {
    log_trace("alpha %.2f, beta %.2f", alpha, beta);
    return pow(x, alpha-1) * pow(1-x, beta-1);
  }

  /**
   * Returns a random variate drawn from the beta distribution
   *
   * Uses the gsl_ran_beta method
   *
   * @param random gsl_rng instance
   */
  static double rvs(const gsl_rng* random) {
    log_trace("alpha %.2f, beta %.2f", alpha, beta);
    return gsl_ran_beta(random, alpha, beta);
  }
};

/**
 *@brief This class implements a Cascade development simulation including
 *       the LPM effect
 *
 * Description of the simulation can be found here: 
 * http://wiki.icecube.wisc.edu/index.php/Cascade_Simulation
 *
 * Cross section formulars used are defined in I3CascadeSimulationCrossSections
 *
 *
 * @author Bernhard Voigt
 */
class I3CascadeSimulation : public I3CascadeDevelopment, public I3CascadeMCCommon::Sampler {

  // define the logger name
  SET_LOGGER("I3CascadeSimulation");

 public:

  /**
   * @brief Constructor
   *
   * @param random random number service from icetray frame
   */
  I3CascadeSimulation(I3RandomServicePtr random);
  
  /**
   * @brief Destructor
   *
   * Frees the gsl objects
   */
  ~I3CascadeSimulation();
   
  virtual void SetRandomNumberGenerator(I3RandomServicePtr);

  /**
   * Simulates a cascade of the given energy
   *
   * @param energy
   */
  void Simulate(I3Particle::ParticleType type, double energy);

  /**
   * Sets the energy threshold down to which particles are tracked (Default is 50 TeV)
   *
   * This shouldn't be lower than 1 GeV, since this is the limit for the Bremsstrahlung 
   * interaction
   *
   * @param thresholdEnergy
   */
  void SetThreshold(double thresholdEnergy) {
    threshold_ = thresholdEnergy;
  }

  /**
   * Returns the threshold energy down to which particles are tracked
   *
   * @return threshold energy
   */
  double GetThreshold() {
    return threshold_;
  }

  /**
   * @brief default energy threshold down to which particles are tracked
   */
  static const double DEFAULT_THRESHOLD;

 private:
  
  I3CascadeSimulation(const I3CascadeSimulation&);

  /**
   * Cacluates the pair productio mean free path for the given energy
   *
   * @param energy
   */
  double PairProductionMeanFreePath(double energy);

  /**
   * Cacluates the bremsstrahlung radiation length (cut off energy is defined as 1 GeV)
   *
   * @param energy
   */
  double BremsstrahlungMeanFreePath(double energy);

  /**
   * Samples the bremsstrahlung cross section to get a fractional energy
   * of a secondary particle produced in an interaction with a incident particle
   * of the given energy
   *
   * @see I3MetropolisHastings::Sample
   *
   * @param energy
   */
  double SampleBremsstrahlungCrossSection(double energy);

  /**
   * Samples the pair production cross section to get a fractional energy
   * of a secondary particle produced in an interaction with a incident particle
   * of the given energy
   *
   * @see I3MetropolisHastings::Sample
   *
   * @param energy
   */
  double SamplePairProductionCrossSection(double energy);

  /**
   * Inits the simulation 
   *
   * Interpolation of the mean free paths is performed
   * The metropolis hastings random samplers are initialized
   */
  void Init();

  /**
   * Inits the metropolis hastings random samplers
   */
  void InitMetropolisHastings();
   
  /**
   * Draw and discard samples from the MH samplers
   * to make them independent of initial conditions
   */
  void BurnIn();

  /**
   * Calcuation and interpolation of the pair production mean free path
   * in an energy range specified by lower and upper (in log10)
   * The interpolation used the given number of points.
   * Default energy range is from 100 GeV to 10 EeV
   *
   * Uses the integration and interpolation routines from gsl
   * 
   * @param lower
   * @param upper
   * @param points
   */
  void InterpolPairProductionMeanFreePath(double lower=2., double upper=13.0, int points=100);

  /**
   * Calcuation and interpolation of the bremsstrahlung radiation length
   * in an energy range specified by lower and upper (in log10)
   * The interpolation used the given number of points.
   * Default energy range is from 100 GeV to 10 EeV
   *
   * Uses the integration and interpolation routines from gsl
   *
   * The integration of the cross section is done down to a lower edge of 
   * I3CascadeSimulation::BREMSSTRAHLUNGCUTOFF
   * this is due to the infrared limit, where the cross section raises to infinite values.
   * 
   * @param lower
   * @param upper
   * @param points
   */
  void InterpolBremsstrahlungMeanFreePath(double lower=2., double upper=13.0, int points=100);

  /**
   * Returns the mean free path of a particle
   *
   * @param particle
   */
  double MeanFreePath(Particle particle);

  /**
   * Propagates a particle
   * 
   * Draws the next interaction point
   * Calls the interaction method to produce secondaries
   *
   * @param particle
   */
  void Propagate(Particle particle);

  /**
   * Creates secondaries produced in an interaction of the given particle
   *
   * @param particle
   */
  void Interaction(Particle particle);

  /**
   * Produces a secondary gamma
   */
  void BremsstrahlungInteraction(Particle electron);

  /**
   * Produces a electron, positron pair
   */
  void PairProductionInteraction(Particle photon);

  /**
   * Deposits the particle
   *
   * Calculates the energy loss profile of the given particle
   * and adds it to the total energy loss profile of this simulation
   */
  void Deposit(Particle particle);

  /**
   * @brief stack holding particles that have to be propagated
   */
  std::stack<Particle> particles_;

  /**
   * @brief counter how many particles have been created
   */
  unsigned int particlesCreated_;

  /**
   * @brief counter how many particles have been deleted
   */
  unsigned int particlesDeleted_;

  /**
   * @brief struct holding beta distribution parameters, provides function 
   * evaluation and random numbers
   */
  beta_dist betaDist_;

  /**
   * @brief random number sampler for low energy bremsstrahlung
   */
  I3MetropolisHastings lowBremsSampler_;

  /**
   * @brief random number sampler for high energy bremsstrahlung
   */
  I3MetropolisHastings highBremsSampler_;

  /**
   * @brief random number sampler for low energy pair production
   */
  I3MetropolisHastings lowPairSampler_;

  /**
   * @brief random number sampler for high energy pair production
   */
  I3MetropolisHastings highPairSampler_;

  /**
   * @brief parametrization of energy loss profile for low energy particles
   */
  I3CascadeParametrization parametrization_;

  /**
   * @brief gsl spline accelerator object used for spline evaluation (bremsstrahlung)
   */
  gsl_interp_accel* interpolAccBrems_;

  /**
   * @brief gsl spline - interpolation of Bresmstrahlung mean free path
   */
  gsl_spline* splineBrems_;

  /**
   * @brief gsl spline accelerator object used for spline evaluation (pair production)
   */
  gsl_interp_accel* interpolAccPair_;

  /**
   * @brief gsl spline - interpolation of Pairproduction mean free path
   */
  gsl_spline* splinePair_;

  /**
   * @brief threshold down to which particles are tracked, default is 50 TeV
   */
  double threshold_;
   
  /**
   * @brief repeat burn-in for each event
   *
   * This ensures that the events are truly independent and can be reproduced
   * given the state of the random number generator.
   */
  bool perEventBurnIn_;

  /**
   * @brief lowest energy of bremsstrahlung photons
   */
  static const double BREMSSTRAHLUNGCUTOFF;

};
#endif






            
