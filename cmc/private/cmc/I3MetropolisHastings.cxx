/**
 * Implementation of Metropolis-Hastings sampler
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3MetropolisHastings.cxx 142695 2016-03-01 18:21:39Z david.schultz $
 *
 * @version $Revision: 142695 $
 * @date $LastChangedDate: 2016-03-01 11:21:39 -0700 (Tue, 01 Mar 2016) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   $LastChangedBy$
 */

// local includes

#include "I3MetropolisHastings.h"
#include "phys-services/I3RandomService.h"
#include "icetray/I3Logging.h"

/* the constructor takes a gsl random instance and functions defined in the beta_dist struct
 * from the I3CascadeSimulation header file. The name parameter is for logging purposes
 */
I3MetropolisHastings::I3MetropolisHastings(I3RandomServicePtr rng,
                                           double (*pdf)(const double&, const double&), 
                                           double (*proposal)(const gsl_rng*),
                                           double (*proposalDistribution)(const double&),
                                           std::string name):
  I3CascadeMCCommon::Sampler(rng),
  pdf_(pdf),
  proposal_(proposal),
  proposalDistribution_(proposalDistribution),
  counter_(0),
  yield_(0),
  name_(name)
{
}

// runs n times the sample method, in order to get rid off the initial state
void I3MetropolisHastings::BurnIn(double energy, double start, unsigned int n) {
  log_debug("Burn in start at %.2f", start);

  // initialize the last state with the given numbers
  x_ = start;
  y_ = pdf_(energy, start);
  propDistY_ = proposalDistribution_(x_);
#ifdef I3_COMPILE_OUT_VERBOSE_LOGGING
  while (n-- > 0) { Sample(energy); }
#else
  double s = -1.0; // initialize to non-sense
  for (; n > 0; n--) {
    // sample steps times
    s = Sample(energy);
  }
  log_debug("Burn in end at %.2f", s);
#endif
}

// sample from pdf(energy), return sample between lower and higher
double I3MetropolisHastings::Sample(double energy, double lower, double higher) {

  // warning counter, if it exceeds a large number, the sampling is canceld
  unsigned int loopcount = 0;
  while (1) {
    counter_++;
    loopcount++;
    if (counter_ % 10000 == 0) {
      // print the efficiency of the sampling, ie. the number of samples yield over the total number of loops
      log_debug("%s - Counter: %lu Efficiency: %.2f", 
               name_.c_str(), counter_, static_cast<double >(yield_)/static_cast<double >(counter_));
    }

    //  draw a candidate from the proposal distribution
    double candidate = proposal_(random_->GSLRng());
    /* if upper and lower bounds are given, the proposal distirbution
     * has to be scaled by scaling the candidate between
     * the upper and lower bounds, this is ok, since it's only a candidate
     * nd the proposal distribution doesn't have to match the target distribution
     * ne could also draw a candidate as long as it's within the bounds, however,
     * his is inefficient. 
     */
    if (candidate < lower || candidate > higher) {
      log_trace("%s - Candidate out of range %.2e", name_.c_str(), candidate);
      candidate = candidate * (higher-lower) + lower;
    }
    log_trace("Candidate %.2f", candidate);

    // calculate proposalDistribtion value, pdf value and acceptance fraction
    double candPropDistY = proposalDistribution_(candidate);
    double probability = pdf_(energy, candidate);
    double acceptance = probability/y_ * propDistY_/candPropDistY;
    log_trace("Acceptance %.2f", acceptance);

    // take the candidate with a probability defined by acceptance
    if (random_->Uniform() < acceptance) {
      // candidate is accepted, keep this state
      x_ = candidate;
      y_ = probability;
      propDistY_ = candPropDistY;
      yield_++; // success -> increase yield count
      return x_;
    }

    // if the sampling is inefficient, warn the user and return the candidate
    // this shouldn't be done, however, if it happens rarely and thus doesn't have an effect on the overall distribution
    if (loopcount >= 1e3) {
      log_warn("%s - Sampling much too inefficient, loopcount limit of 1e3 reached! Energy %.2f", name_.c_str(), energy);
      return candidate;
    }
  }
}


