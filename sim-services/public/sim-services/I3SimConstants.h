/**
 * copyright  (C) 2012
 * the icecube collaboration
 * @version $Id: I3SimConstants.h $
 * @file I3SimConstants.h
 * @date $Date: $
 */

#ifndef I3SIMCONSTANTS_H_INCLUDE
#define I3SIMCONSTANTS_H_INCLUDE

#include <icetray/I3Units.h>
#include <dataclasses/I3Constants.h>
#include <dataclasses/physics/I3Particle.h>
#include <cmath>

/**
 * @brief A list of static variables commonly used by physics simulation
 *        specifically cross-section calculations
 *
 * Just a namespace filled with constants -- add a line
 * <tt> using namespace I3SimConstants </tt> to get access to them
 * directly, or just use @c I3SimConstants::pi, for example.
 */
namespace I3SimConstants
{
  /**
   *  Particle Masses, etc. taken from
   *  http://pdg.lbl.gov/2012/reviews/rpp2012-rev-phys-constants.pdf
   *  unless otherwise stated
   */
   
  /**
   *  Mass of electron
   */
  
  static const double m_e = 0.510998928 * I3Units::MeV / (I3Constants::c * I3Constants::c);
  
  /**
   *  Mass of muon
   *  http://pdg.lbl.gov/2012/listings/rpp2012-list-muon.pdf
   */ 
  
  static const double m_mu = 105.6583715 * I3Units::MeV / (I3Constants::c * I3Constants::c);
   
  /**
   *  Mass of tau
   *  http://pdg.lbl.gov/2012/listings/rpp2012-list-tau.pdf
   */ 
  
  static const double m_tau = 1776.82 * I3Units::MeV / (I3Constants::c * I3Constants::c);
  
  /**
   *  Mass of proton
   */
  
  static const double m_p = 938.272046 * I3Units::MeV / (I3Constants::c * I3Constants::c);
  
  /**
   * Mass of neutron
   */
   
  static const double m_n = 939.565379 * I3Units::MeV / (I3Constants::c * I3Constants::c);
   
  /**
   *  Mass of W boson
   */
  
  static const double m_W = 80.385 * I3Units::GeV / (I3Constants::c * I3Constants::c);
   
  /**
   *  Mass of Z boson
   */
  
  static const double m_Z = 91.1876 * I3Units::GeV / (I3Constants::c * I3Constants::c);

  /**
   *  Planck Constant
   */
  
  static const double h = 6.62606957e-34 * I3Units::joule * (I3Units::second);

  /**
   *  reduced Planck Constant
   */

  static const double hbar = h / (2. * I3Constants::pi);

  /**
   *  unit conversion constant (hbar * c)
   */

  static const double hbarc = hbar * I3Constants::c;
     
  /**
   *  Fermi coupling constant
   *  / (hbar * c)^3 is assumed
   */

  static const double G_Fermi = 1.1663787e-5 / (I3Units::GeV * I3Units::GeV);
  
  /**
   *  Permittivity of free space
   */
  
  static const double epsilon_0 = 8.854187817e-12 * I3Units::coulomb / (I3Units::V * I3Units::m);

  /**
   *  Permeability of free space
   */
  
  static const double mu_0 = 4 * I3Constants::pi * 10e-7 * I3Units::V * I3Units::second / (I3Units::A * I3Units::A);
   
  /**
   *  Fine-structure constant
   */

  static const double alpha = 1./137.03599976; // (I3Units::eSI * I3Units::eSI) / ( 4 * I3Constants::pi * epsilon_0 * hbarc);

  /**
   *  Boltzmann constant
   */
     
  static const double k_B = 1.3806488 * I3Units::joule / I3Units::kelvin;
   
  /**
   *  Weak-mixing Angle
   *  \f$sin^{2} \hat{\theta} (M_{Z})
   *  cos \theta_{W}\f$
   */
     
  static const double sinsq_theta_W = 0.23116;
   
  static const double cos_theta_W = m_W / m_Z;
  
  
  /**
   *  Cabibbo Angle
   */
     
  static const double sin_theta_C = 0.4561;
   
  static const double cos_theta_C = 0.9746;
  
  /** 
   *  Shower development parameters, taken from IceCube IR icecube/201210001
   */
  struct ShowerParameters {
      /** 
       * @param[in] type    Type of particle that initiates the shower
       * @param[in] energy  Energy of the primary particle
       * @param[in] density Density of the medium
       */
      ShowerParameters(I3Particle::ParticleType type, double energy,
          double density=0.9216*(I3Units::g/I3Units::cm3));
      /** Shape parameter of a gamma distribution (dimensionless) */
      double a; 
      /** Scale length of a gamma distribution (meters) */
      double b;
      /** Ratio of total cherenkov track yield relative to EM cascade */
      double emScale;
      /** Standard deviation of relative flucuations in the shower light yield */
      double emScaleSigma;
  };

};

#endif //I3SIMCONSTANTS_H_INCLUDED

