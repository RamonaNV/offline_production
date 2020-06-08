/**
 * This file is the header file of the
 * I3CascadeMCModule class.
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeMCModule.h 43507 2008-03-20 09:32:11Z bvoigt $
 *
 * @version $Revision: 43507 $
 * @date $LastChangedDate: 2008-03-20 05:32:11 -0400 (Thu, 20 Mar 2008) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: bvoigt $
 */

#ifndef I3_CASCADE_MC_SERVICE_H_INCLUDED
#define I3_CASCADE_MC_SERVICE_H_INCLUDED

class I3CascadeSplit;
class I3CascadeMuonSplit;

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "phys-services/I3RandomService.h"
#include "sim-services/I3PropagatorService.h"

/**
 * @brief Cascade Monte Carlo Service
 * 
 * This module provides the functionality
 * to replace cascade-like in the MCTree with
 * with a list of sub-cascades to simulate
 * the cascade longitudinal development.
 */
class I3CascadeMCService : public I3PropagatorService {

  SET_LOGGER("I3CascadeMCService");

 public:

  ~I3CascadeMCService() {}

  I3CascadeMCService(I3RandomServicePtr r);
  
  virtual std::vector<I3Particle> Propagate(I3Particle& p, DiagnosticMapPtr frame, I3FramePtr);
  virtual void SetRandomNumberGenerator(I3RandomServicePtr random);

  void Simulate(I3Particle& c, std::vector<I3Particle>& d);

  void SetEnergyThresholdMuons( double t);

  void SetMaxMuons( int i );

  void SetThresholdSplit(double t );

  void SetEnergyThresholdSimulation(double t);
  
  void SetStepWidth(int w );
 
 private:

 /**
   * @brief random service provided by the frame
   */
  I3RandomServicePtr random_;

  /**
   * @brief hadronic particles with an energy exceeding this 
   * threshold will generate muons.
   */
  double energyThresholdMuons_;

  /**
   * @brief maximal number of generated muons.
   */
  int maxMuons_;

  /**
   * @brief particles with an energy exceeding this threshold will be split
   */
  double energyThresholdSplit_;

  /**
   * @brief particles with an energy exceeding this threshold will be split 
   * and the the longitudinal development is simulated
   */
  double energyThresholdSimulation_;

  /**
   * @brief distance between individual sub-cascades in units of radiation length
   */
  int stepWidth_;

   /**
   * @brief cascade splitting class
   */
  boost::shared_ptr<I3CascadeSplit> splitter_;

  /**
   * @brief cascade moun splitting class
   */
  boost::shared_ptr<I3CascadeMuonSplit> muonSplitter_;

};

#endif
