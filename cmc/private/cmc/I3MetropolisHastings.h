/**
 * A metropolis hasting sampler for bremsstrahlung and pair production cross sections 
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3MetropolisHastings.h 143040 2016-03-11 17:35:25Z nega $
 *
 * @version $Revision: 143040 $
 * @date $LastChangedDate: 2016-03-11 10:35:25 -0700 (Fri, 11 Mar 2016) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   $LastChangedBy$
 */

#ifndef I3METROPOLISHASTINGS_h_
#define I3METROPOLISHASTINGS_h_

// icetray includes
#include "icetray/I3Logging.h"
#include "cmc/I3CascadeMCCommon.h"

// gsl includes
#include <gsl/gsl_rng.h>

// c++ includes
#include <string>

/**
 * @brief This class implements a Metropolis Hastings sampler
 *
 * It uses a pdf and a proposal distribution to efficiently
 * produce random samples from the pdf
 *
 * The proposal distribution used is a beta distribution, which
 * matches almost the shape of the bremsstrahlung and pair production
 * cross sections, using the right parameters. See I3CascadeSimulation
 * for the usage.
 *
 * For details on the algorithm, take a look at :
 * http://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm
 */
class I3MetropolisHastings : public I3CascadeMCCommon::Sampler {

  // set logging name for icetray
  SET_LOGGER("I3MetropolisHastings");

 public:

  /**
   * Default Constructor
   */
  I3MetropolisHastings() {}

  /**
   * Constructor
   *
   * @param rng      a gsl_random instance
   * @param pdf      a function with two parameters, ie. the energy and
   *                 fractional energy of the cross section
   * @param proposal a gsl random instance sampling from the proposal
   *                 distribution which is a beta distribution
   * @param proposalDistribution a beta function in x
   * @param name     name of this instance for logging purposes
   *
   * The proposal gsl random instance as well as the proposalDistribution are implemented 
   * in a struct beta_dist in the header file I3CascadeSimulation.h
   */
  I3MetropolisHastings(I3RandomServicePtr rng,
                       double (*pdf)(const double&, const double&), 
                       double (*proposal)(const gsl_rng*),
                       double (*proposalDistribution)(const double&),
                       std::string name);

  /**
   * Get the total sample loop count 
   */
  unsigned long GetCounter() {
    return counter_;
  }

  /**
   * Reset the total sample loop count
   */
  void ResetCounter() {
    counter_ = 0;
  }

  /**
   * Get the number of returned samples
   */
  unsigned long GetYield() {
    return yield_;
  }

  /**
   * Reset the number of returned samples
   */
  void ResetYield() {
    yield_ = 0;
  }


  /**
   * Samples n times starting at start position
   *
   * Used to forget the initial state
   *
   * @param energy parameter of underlying pdf
   * @param start starting position
   * @param n number of samples
   */
  void BurnIn(double energy, double start, unsigned int n=100);

  /**
   * Samples the underlying pdf in the given range
   *
   * @param energy parameter to the underlying pdf
   * @param lower edge of sample range
   * @param higher edge of sample range
   * @return random number
   */
  double Sample(double energy, double lower=0, double higher=1);

 private:

  /**
   * @brief pdf, ie a cross section with energy and fractional energy parameters
   */
  double (*pdf_)(const double&, const double&);

  /**
   * @brief a gls random instance drawing samples from the proposal distribution
   */
  double (*proposal_)(const gsl_rng*);

  /**
   * @brief a proposal distribution, in this case a beta distribution in x
   */
  double (*proposalDistribution_)(const double&);

  /**
   * @brief loop counter
   */
  unsigned long counter_;

  /**
   * @brief counts how many successful attemps have been made to draw a new random number
   */
  unsigned long yield_;

  /**
   * @brief name of this instance, for logging purpose
   */
  std::string name_;

  /**
   * @brief last drawn sample
   */
  double x_;

  /**
   * @brief last value of pdf(x)
   */
  double y_;

  /**
   * @brief last value of proposalDist(x)
   */
  double propDistY_;

};

#endif



