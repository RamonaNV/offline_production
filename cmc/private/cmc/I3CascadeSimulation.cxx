/**
 *@brief Implementation I3CascadeSimulation class
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSimulation.cxx 150508 2016-10-04 14:44:25Z jvansanten $
 *
 * @version $Revision: 150508 $
 * @date $LastChangedDate: 2016-10-04 08:44:25 -0600 (Tue, 04 Oct 2016) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: jvansanten $
 */

#include <algorithm>
#include <cmath>
#include <iostream>

// local includes
#include "cmc/I3CascadeSimulation.h"
#include "cmc/I3CascadeDevelopment.h"
#include "cmc/I3CascadeParametrization.h"
#include "cmc/I3CascadeSimulationCrossSections.h"
#include "cmc/I3MetropolisHastings.h"

// Gnu Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_histogram.h>

// limits for numbers
#include <climits>

// IceTray includes
#include "phys-services/I3RandomService.h"
#include "phys-services/I3GSLRandomService.h"
#include "icetray/I3Units.h"

/**
 * Local definition of a wrapper function needed to call the gsl_integration routine
 * 
 * @param params energy fraction
 * @param y      double value
 */
double BremsstrahlungCrossSectionWrapper(double y, void* params) {
  double energy = *(double*) params;
  return I3CascadeSimulationCrossSections::BremsstrahlungCrossSection(energy, y);
}

/**
 * Local definition of a wrapper function needed to call the gsl_integration routine
 * 
 * @param params energy fraction
 * @param y      double value
 */
double PairProductionCrossSectionWrapper(double y, void* params) {
  double energy = *(double*) params;
  double result  = I3CascadeSimulationCrossSections::PairProductionCrossSection(energy, y);
  //log_debug("eval pp(%.3e, %.3f: %.3e", energy, y, result);
  return result;
}

// definition of class constants
const double I3CascadeSimulation::DEFAULT_THRESHOLD = 50 * I3Units::TeV;
const double I3CascadeSimulation::BREMSSTRAHLUNGCUTOFF = 1 * I3Units::GeV;

// public fields of struct must be defined
double beta_dist::alpha;
double beta_dist::beta;

/**
 * Constructor with random service given
 *
 * same as above, except random is not created
 */
I3CascadeSimulation::I3CascadeSimulation(I3RandomServicePtr random) :
  I3CascadeDevelopment(),
  I3CascadeMCCommon::Sampler(random),
  interpolAccBrems_(0),
  splineBrems_(0),
  interpolAccPair_(0),
  splinePair_(0),
  threshold_(I3CascadeSimulation::DEFAULT_THRESHOLD),
  perEventBurnIn_(false)
{ 
  Init();
}

void
I3CascadeSimulation::SetRandomNumberGenerator(I3RandomServicePtr rng)
{
  random_ = rng;
  lowBremsSampler_.SetRandomNumberGenerator(rng);
  highBremsSampler_.SetRandomNumberGenerator(rng);
  lowPairSampler_.SetRandomNumberGenerator(rng);
  highPairSampler_.SetRandomNumberGenerator(rng);
  // If someone set a random number service, they want complete
  // control over event-to-event correlations in state.
  perEventBurnIn_ = true;
}

/**
 * Destructor - frees all gsl struct that have been initialized
 */
I3CascadeSimulation::~I3CascadeSimulation() {
  if (interpolAccBrems_ != 0) {
    gsl_interp_accel_free(interpolAccBrems_);
  }
  if (splineBrems_ != 0) {
    gsl_spline_free(splineBrems_);  
  }
  if (interpolAccPair_ != 0) {
    gsl_interp_accel_free (interpolAccPair_);
  }
  if (splinePair_ != 0) {
    gsl_spline_free(splinePair_);
  }
}

/**
 * Creates different metropolis hasting samplers for pair production
 * and bremsstrahlung cross sections
 *
 * The parameters of the proposal distribution (beta distribution) are
 * set individualy in the sample methods below.
 *
 * There are different instances of the sampler since the sampler keeps
 * the old state and they shouldn't mix accross the different cross sections
 * for high and low energies, as well as for the different proccesses.
 */
void I3CascadeSimulation::InitMetropolisHastings() {
  log_trace("Initialization of MetropolisHastings sampler");

  // creata a beta distribution instance
  beta_dist betaDist_ = beta_dist();
  // set parameter of beta distribution, dummy values ony for the burn in
  betaDist_.alpha = .5;
  betaDist_.beta = .5;

  // metropolis hastings sampler for bremsstrahlung for low and high energy
  lowBremsSampler_ = I3MetropolisHastings(random_,
                                          I3CascadeSimulationCrossSections::BremsstrahlungCrossSection,
                                          betaDist_.rvs,
                                          betaDist_.pdf, "lowBrems");
  highBremsSampler_ = I3MetropolisHastings(random_,
                                           I3CascadeSimulationCrossSections::BremsstrahlungCrossSection,
                                           betaDist_.rvs,
                                           betaDist_.pdf, "highBrems");

  // metropolis hastings sampler for pairproduction forlow and high energy
  lowPairSampler_ = I3MetropolisHastings(random_,
                                         I3CascadeSimulationCrossSections::PairProductionCrossSection,
                                         betaDist_.rvs,
                                         betaDist_.pdf, "lowPair");
  highPairSampler_ = I3MetropolisHastings(random_,
                                          I3CascadeSimulationCrossSections::PairProductionCrossSection,
                                          betaDist_.rvs,
                                          betaDist_.pdf, "highPair");
  BurnIn();
}

void I3CascadeSimulation::BurnIn()
{
  betaDist_.alpha = .5;
  betaDist_.beta = .5;
  lowBremsSampler_.BurnIn(1 * I3Units::TeV, .1, 100);
  highBremsSampler_.BurnIn(1 * I3Units::PeV, .1, 100);
  lowPairSampler_.BurnIn(1 * I3Units::TeV, .1, 100);
  highPairSampler_.BurnIn(1 * I3Units::PeV, .1, 100);
}

/**
 * integrate and interpolate pair production cross section to get the mean free path
 */
void I3CascadeSimulation::InterpolPairProductionMeanFreePath(double lower, double upper, int points) {
  log_debug("Interpolation of total Pairproduction crosssection");

  // x, y points
  double x[points];
  double y[points];
  
  // integration boundaries
  double a = 0.0;
  double b = 1.0;
  // integration result and error
  double result, error;

  // function to integrate, parameter to function is set in the loop
  gsl_function fct;
  fct.function = &PairProductionCrossSectionWrapper;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);

  // calculate the total cross section
  for(int i=0; i < points; i++) {
    // evenly distributed points in powers of GeV 
    double energy = pow(10, lower + i * (upper-lower)/(points-1));
    // set function parameter
    fct.params = &energy;
    // integrate using gsl routine, docs here:
    // http://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
    gsl_integration_qags(&fct, a, b, 0, 1e-5, 1000, w, &result, &error);
    double yi = result;
    //log_debug("total pp cs: %.3e", result);
    x[i] = energy;
    y[i] = 1/yi; // mean free path is 1/totalCrossSection
  }
  gsl_integration_workspace_free(w);

  // interpolate the mean free path
  interpolAccPair_ = gsl_interp_accel_alloc();
  splinePair_ = gsl_spline_alloc(gsl_interp_cspline, points);
  gsl_spline_init(splinePair_, x, y, points);
}

/**
 * integrate and interpolate bremsstrahlung cross section to get mean free path
 */
void I3CascadeSimulation::InterpolBremsstrahlungMeanFreePath(double lower, double upper, int points) {
  log_debug("Interpolation of total Bremsstrahlung crosssection");

  // x, y points
  double x[points];
  double y[points];
  
  // integration boundaries
  double a; // individually set for each energy
  double b = 1.0;
  // integration result and error
  double result, error;

  // function to integrate, params are set in the loop
  gsl_function fct;
  fct.function = &BremsstrahlungCrossSectionWrapper;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  // if lower end of energy region is not larger than 1% of the bremsstrahlung cutoff, force it to be so
  // this is to omit the infrared singularity in the integration
  if (BREMSSTRAHLUNGCUTOFF/pow(10,lower) > .01) {
    lower = log10(BREMSSTRAHLUNGCUTOFF * 100);
    log_warn("Bremsstrahlung cutoff too large for lower edge of interpolation region! Set lower energy to energy %.2f", pow(10,lower));
  }

  // calculate the total cross section
  for(int i=0; i < points; i++) {
    // evenly distributed points in log 
    double energy = pow(10, lower + i * (upper-lower)/(points-1));
    // set function parameter
    fct.params = &energy;
    // lower energy fraction is the energy fraction of the bremsstrahlung cut off of the energy of the particle
    a = BREMSSTRAHLUNGCUTOFF/energy;
    // integrate using gsl routine, docs here: 
    // http://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
    // numbers are the rel and abs error, respectively
    gsl_integration_qags(&fct, a, b, 1e-5, 1e-5, 1000, w, &result, &error);
    double yi = result;
    //log_debug("total bs cs: %.3e", result);
    x[i] = energy;
    y[i] = 1/yi; // mean free path is 1/totalCrossSection
  }
  gsl_integration_workspace_free(w);

  // interpolate the mean free path
  interpolAccBrems_ = gsl_interp_accel_alloc();
  splineBrems_ = gsl_spline_alloc(gsl_interp_cspline, points);
  gsl_spline_init(splineBrems_, x, y, points);
}

/**
 * Calls the interpolation routine and the mh-sampler initialization
 */
inline void I3CascadeSimulation::Init() {
  log_trace("Initialization of simulation object");

  // inits all spline and interpolation pointers
  InterpolPairProductionMeanFreePath();
  InterpolBremsstrahlungMeanFreePath();

  // set threshold for parametrization of energy loss profile of low energy particles
  parametrization_.SetThreshold(1 * I3Units::GeV);

  // create and burn in off metropolis-hastings samplers
  InitMetropolisHastings();
}

// return pair productio mean free path from interpolation
double I3CascadeSimulation::PairProductionMeanFreePath(double energy) {
  double res =  gsl_spline_eval(splinePair_, energy, interpolAccBrems_);
  log_trace("SplPair_(%.3e): %.3e\n", energy, res);
  return res;
}
// return bremsstrahlung mean free path from interpolation
double I3CascadeSimulation::BremsstrahlungMeanFreePath(double energy) {
  double res = gsl_spline_eval(splineBrems_, energy, interpolAccBrems_);
  log_trace("SplBrems_(%.3e): %.3e\n", energy, res);
  return res;
}

// sample from the bremsstrahlung cross section
double I3CascadeSimulation::SampleBremsstrahlungCrossSection(double energy) {
  log_trace("Sampling Bremsstrahlung");

  // set parameter to beta fct depending on energy
  // this is to tune the proposal distribution, that it matches the shape
  // of the cross section, the numbers used are obtained from fits of the beta
  // distribution to the differential cross section.
  if (energy < 1 * I3Units::PeV) {
    betaDist_.alpha = .35;
    betaDist_.beta = 1.5;
    // sample with a lower edge given by the fraction of bremsstrahlung cut-off and the particles energy
    return lowBremsSampler_.Sample(energy, BREMSSTRAHLUNGCUTOFF/energy);
  } else {
    betaDist_.alpha = .45;
    betaDist_.beta = .75;
    return highBremsSampler_.Sample(energy, BREMSSTRAHLUNGCUTOFF/energy);
  }
}

// sample from the pair production cross section
double I3CascadeSimulation::SamplePairProductionCrossSection(double energy) {
  log_trace("Sampling Pairproduction");
  // set parameter to beta fct depending on energy
  // this is to tune the proposal distribution, that it matches the shape
  // of the cross section, the numbers used are obtained from fits of the beta
  // distribution to the differential cross section.
  if (energy < 5 * I3Units::PeV) {
    betaDist_.alpha = .85;
    betaDist_.beta = .85;
    return lowPairSampler_.Sample(energy);
  } else {
    betaDist_.alpha = .4;
    betaDist_.beta = .4;
    return highPairSampler_.Sample(energy);
  }
}

// run the simulation
// a first particle is injected, brought to interaction, secondaries are
// put into the particle list, if the energy of secondaries is still above
// the threshold, run the simulation for these
void I3CascadeSimulation::Simulate(I3Particle::ParticleType type, double energy) {
  log_trace("Entering Simulate");

  // ensure that the samples depend only the state of the RNG at
  // the point of call, and not the order of events previously simulated
  if (perEventBurnIn_)
    BurnIn();

  // check particle stack
  if (particles_.size())
    log_fatal("%zu particles are still on the propagation stack from the last run.", particles_.size());
  // init fresh energy loss profile, that's the total energy loss profile of this shower
  energyLossProfile_.clear();
  energyLossProfile_.reserve(200);

  // incident particle is an electron an position 0
  particles_.push(Particle(energy, 0.0, Particle::ELECTRON));
  
  // propaget all particles to track until the list is empty
  long counter = 0; // loop counter for debugging
  while(particles_.size() > 0) {
    if (counter % 1000 == 0) {
      log_debug("Simulation loop count: %li", counter);
    }
    // get the last particle and delete it from the vector
    // the object will be deleted in the Propaget/Interaction methods
    Particle particle = particles_.top();
    particles_.pop();

    // propaget the particle
    Propagate(particle);
    counter++;
  }
}

// get the particles interaction mean free path
double I3CascadeSimulation::MeanFreePath(Particle particle) {
  return (particle.type == Particle::ELECTRON) ? BremsstrahlungMeanFreePath(particle.energy)
    : PairProductionMeanFreePath(particle.energy);
}

// draw interaction point, call interaction routine
void I3CascadeSimulation::Propagate(Particle particle) {
  log_trace("Entering Propagate");
  double lambda = MeanFreePath(particle) * I3Units::cm;
  // random number acoording to the interaction probability
  // P(x) = 1-exp(x/lambda)
  double x = -log(random_->Uniform()) * lambda;
  // particle reached the new point
  particle.x += x;
  //simulate interaction
  Interaction(particle);
}

// calls brems or pair production interaction
void I3CascadeSimulation::Interaction(Particle particle) {
  log_trace("Entering Interaction");
  if (particle.type == Particle::ELECTRON)
    return BremsstrahlungInteraction(particle);
  else
    return PairProductionInteraction(particle);
}

// bremsstrahlung interactin, creates photon, deposits particles when energy is below threshold
void I3CascadeSimulation::BremsstrahlungInteraction(Particle electron) {
  // draw fractional energy from cross section
  double energyFraction = SampleBremsstrahlungCrossSection(electron.energy);
  // create photon with fractional energy at this interaction point
  // the memory allocated by new is released in Deposit which is called for every particle
  Particle photon(electron.energy * energyFraction, electron.x, Particle::PHOTON);

  // substract photon energy from electron
  electron.energy -= photon.energy;

  // particle will be tracked only if their energies exceed the threshold energy
  if (electron.energy > threshold_) {
    particles_.push(electron);
  } else {
    // deposit particle according to the energy deposit parametrization and delete it
    Deposit(electron);
  }

  // photon will be tracked if energy above threshold
  if (photon.energy > threshold_) {
    particles_.push(photon);
  } else {
    // deposit particle according to the energy deposit parametrization and delete it
    Deposit(photon);
  }
}

// pair production interaction, creates electron and positron, deletes photon, deposits particles when
// their energy is below threshold
void I3CascadeSimulation::PairProductionInteraction(Particle photon) {
  // draw fractional energy from cross section
  double energyFraction = SamplePairProductionCrossSection(photon.energy);
  // create electron pair with total photon energy at the interaction point
  Particle electron1(photon.energy * energyFraction, photon.x, Particle::ELECTRON);
  Particle electron2(photon.energy - electron1.energy, photon.x, Particle::ELECTRON);

  // particles will be tracked only if their energies exceeds the threshold energy
  if (electron1.energy > threshold_) {
    particles_.push(electron1);
  } else {
    // deposit particle according to the energy deposit parametrization and delete it
    Deposit(electron1);
  }

  // particles will be tracked only if their energies exceeds the threshold energy
  if (electron2.energy > threshold_) {
    particles_.push(electron2);
  } else {
    // deposit particle according to the energy deposit parametrization and delete it
    Deposit(electron2);
  }
}

// calcuates energy deposit profile for low energy particles
// adds individual profile to the total profile
// why isn't this a boost::shared_ptr or reference?
void I3CascadeSimulation::Deposit(Particle particle) {
  if( ! std::isnormal(particle.energy) ){
    log_error("Particle energy is not normal %e", particle.energy);
    return;
  }
  // get the position, cast to radiation length (binning is in radiation length)
  int depth = static_cast<int >(floor(particle.x / RADIATION_LENGTH));

  // if particle energy is larger than parametrization threshold
  if (particle.energy > parametrization_.GetThreshold()) {
    parametrization_.Simulate(I3Particle::EMinus, particle.energy);

    // get length of the shower
    int length = static_cast<int >(parametrization_.GetLength());

    // resize energy loss profile vector, if it can't take the new profile
    if (static_cast<int >(energyLossProfile_.size()) < depth + length) {
      energyLossProfile_.resize(depth + length, 0.0);
    }

    // loop until the end of the shower and get energy deposit
    for(int i = 0; i < length; i++) {
      // get the energy loss between this step and the next
      double energy = parametrization_.GetEnergyDeposit(i, i+1);
      // add energy loss to the total energy loss of the shower
      energyLossProfile_[depth+i] += energy;
    }

  } else {
    // just add the particle's energy at the depth, increase length if neccesary
    if (static_cast<int >(energyLossProfile_.size()) < depth + 1) {
      energyLossProfile_.resize(depth + 1, 0.0);
    }
    energyLossProfile_[depth] += particle.energy;
  }
}
