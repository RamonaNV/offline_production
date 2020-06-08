/**
 * Bremsstrahlung and Pairproduction cross sections including LPM suppresion
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSimulationCrossSections.h 48260 2008-08-18 13:10:15Z bvoigt $
 *
 * @version $Revision: 48260 $
 * @date $LastChangedDate: 2008-08-18 07:10:15 -0600 (Mon, 18 Aug 2008) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: bvoigt $
 */

#ifndef I3_CASCADE_SIMULATION_CROSSSECTIONS_H_INCLUDED
#define I3_CASCADE_SIMULATION_CROSSSECTIONS_H_INCLUDED

// IceTray includes
#include "icetray/I3Units.h"
#include "icetray/I3Logging.h"

/**
 *@brief This class implements Bremsstrahlung and Pairproduction cross sections including 
 *       suppresion due to the LPM effect
 *
 * The cross section are take from Spencer Klein's review:
 * "Suppression of bremsstrahlung and pair production due to environmental factors"
 * http://arxiv.org/pdf/hep-ph/9802442
 *
 * @author Bernhard Voigt
 */
class I3CascadeSimulationCrossSections
{
  
  // define the logger name
  SET_LOGGER("I3CascadeSimulationCrossSections");
  
 public:
  /**
   * Caclulates the bremsstrahlung differential cross section for a given energy
   * and given energy fraction of the secondary particle
   */
  static double BremsstrahlungCrossSection(const double& energy, const double& y);

  /**
   * Caclulates the pair production differential cross section for a given energy
   * and given energy fraction of the secondary particle
   */
  static double PairProductionCrossSection(const double& energy, const double& y);

 private:

  /**
   * Calculates the mean Z of a component of different elements
   * 
   * @param exp is the power the sum of components can be raised before
   * the mean is calcualted (a+b+c)**power/3
   */
  static double MeanZ(const double& exp);

  /**
   * Function used in the cross section calculation
   */
  static double Xi(const double& s);

  /**
   * Function used in the cross section calculation
   */
  static double Phi(const double& s);

  /**
   * Function used in the cross section calculation
   */
  static double Psi(const double& s);

  /**
   * Function used in the cross section calculation
   */
  static double G(const double& s);


  /**
   * @brief number of elements of the material
   */
  static const int NUMBER_OF_COMPONENTS;

  /**
   * @brief atomic charge of material, filled in the constructor with values
   *        for wather this is (1, 1, 8) (H_2 O), set in the implemenatin file
   */
  static const int Z[];

  /**
   * @brief radiation length in material (here water/ice)
   */
  static const double X0; // Unit is cm

  /**
   * @brief fine structur constant
   */
  static const double ALPHA;

  /**
   * @brief electron radius
   */
  static const double ELECTRON_RADIUS; // Unit is cm

  /**
   * @brief Avogadro's number
   */
  static const double N; // Unit is 1/mol

  /**
   * @brief Density of material (here ice)
   */
  static const double DENSITY; // Unit g/cm**3

  /**
   * @brief atomic weight of water molecule
   */
  static const double A; // Unit g/mol

  /**
   * @brief ratio of Z/A, the atomic charge over atomic weight (here for water)
   */
  static const double Z_A_RATIO;

  /**
   * @brief energy threshold for which the LPM suppression sets in
   */
  static const double ENERGY_THRESHOLD; // this is 7.7e3 GeV * X0 [g/cm^2]

  /**
   * @brief number of nucleons in the material seen by an electron/photon in the 
   *        interaction. Taken from the reference paper and particle data booklet
   *        The value is set in the implementation file
   */
  static const double NUMBER_OF_NUKLEONS;
  
};
#endif
