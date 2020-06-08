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

#ifndef I3_CASCADE_MC_COMMON_H_INCLUDED
#define I3_CASCADE_MC_COMMON_H_INCLUDED

#include "icetray/I3PointerTypedefs.h"
	 
#include "dataclasses/physics/I3Particle.h"
I3_FORWARD_DECLARATION(I3RandomService);

namespace I3CascadeMCCommon {

  /**
   * @brief checks wether the particle is cascade like
   * 
   * @param particle I3Particle
   */
  bool IsCascade(const I3Particle& particle);

  /**
   * @brief checks wether the particle is hadron shower
   * 
   * @param particle I3Particle
   */
  bool IsHadron(const I3Particle& particle);

  /**
   * @brief Scales hadronic shower energy
   *
   * At high energies hadronic showers are becoming more
   * electromagnetic like and hence the constant scaling with .8
   * applied in Photonics is not suitable.
   * Here a energy dependent scaling factor is applied
   * Taken from M. Kowalski's internal note on the simulation of Hadr. Cascades
   *
   * @param energy energy of the hadron shower
   * @param rand
   * @return rescaled energy
   */ 
  double HadronEnergyCorrection(double energy, boost::shared_ptr<I3RandomService> rand);
  
  /**
   * @brief a mix-in class for class that hold a random number generator
   */
  class Sampler {
  public:
    Sampler(I3RandomServicePtr rng=I3RandomServicePtr()) : random_(rng) {}
    virtual void SetRandomNumberGenerator(I3RandomServicePtr rng) { random_ = rng; }
  protected:
    I3RandomServicePtr random_;
  };
 
};

#endif
