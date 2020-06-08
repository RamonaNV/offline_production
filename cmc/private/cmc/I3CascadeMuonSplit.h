/**
 * Generation of muons for a realistic treatment 
 * of the longitudinal development of cascasdes
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSplit.h 33096 2007-06-11 15:09:40Z bvoigt $
 *
 * @version $Revision: 33096 $
 * @date $LastChangedDate: 2007-06-11 17:09:40 +0200 (Mon, 11 Jun 2007) $
 * @author Sebastian Panknin <panknin@physik.hu-berlin.de>   
 * Last changed by: $LastChangedBy: bvoigt $
 */

#ifndef I3CASCADEMUONSPLIT_h_
#define I3CASCADEMUONSPLIT_h_

//Standard C/C++ includes
#include <cmath>
#include <vector>

// Local includes
// #include "cmc/I3CascadeDevelopment.h"
// #include "cmc/I3CascadeParametrization.h"
// #include "cmc/I3CascadeSimulation.h"
#include "cmc/I3CascadeMCCommon.h"

// icetray includes
#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"

/**
 *@brief This class defines methods to generate muons from a cascade
 *
 * There is only one method of interest to the user which is the 
 * split method.
 * It takes a I3Particle of a cascade typ, generates a set of muons,
 * reduces the energy of the cascade and returns the corresponding 
 * I3Particles for the muons.
 *
 * @author Sebastian Panknin
 */
class I3CascadeMuonSplit : public I3CascadeMCCommon::Sampler {

  // set logging name for icetray
  SET_LOGGER("I3CascadeMuonSplit");

 public:

  /**
   * Creates a I3CascadeMuonSplit object
   *
   * @param random service pointer
   */
  I3CascadeMuonSplit(I3RandomServicePtr random);

  /**
   * Default destructor
   * 
   */
  virtual ~I3CascadeMuonSplit();
  
  /**
   * @brief Generate Muons for a given I3Particle
   *
   * For a (hadronic) cascade-like particle muons according to
   * the cascade energy are generated. The number and energies 
   * of the muons is taken from a parameterization.
   * The muons have the same origian and the same direction as
   * the orgriginal cascade. The are positivly charged 
   * (I3Particle::MuPlus).
   * The energy of the cascade is reduced by the energy 
   * of the muons.
   *
   * @param cascade I3Particle
   * @return vector of I3Particles 
   */
  std::vector<I3Particle> GenerateMuons(I3Particle &cascade);

  /**
   * @brief set the energyThresholdMuons
   */
  void SetEnergyThresholdMuons(double threshold);
  /**
   * @brief set the maxMuons number
   */
  void SetMaxMuons(int max);

  static const double DEFAULT_ENERGY_THRESHOLD;
  static const int DEFAULT_MAX_MUONS;

 private:
  I3CascadeMuonSplit(const I3CascadeMuonSplit& cascadeMuonSplit);
  
  /**
   * @brief energy treshold for muons to consider
   */
  double energyThresholdMuons_;

  /**
   * @brief maximal number of muons to consider
   */
  int maxMuons_;


};
#endif
