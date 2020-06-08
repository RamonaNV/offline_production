#include <I3Test.h>

#include <vector>
#include <limits>
#include "phys-services/I3SPRNGRandomService.h"
#include "phys-services/I3SPRNGRandomServiceFactory.h"
#include "dataclasses/physics/I3Particle.h"
#include "icetray/I3Units.h" 
#include "cmc/I3CascadeSplit.h"


// checks the direction, position and total energy of sub-cascade list
void Ensure(I3Particle p, std::vector<I3Particle> subCascades) {

  std::vector<I3Particle>::const_iterator iter = subCascades.begin();

  // check location, this is a float comparison
  ENSURE_DISTANCE((iter->GetPos() - p.GetPos()).Magnitude(), 0.5*0.358 / 0.9216, 1e-7,
      "Start location of sub-cascades is 1/2 radiation length from original cascade");

  // check direction and sum energies
  double totalEnergy = 0.0;
  for(; iter != subCascades.end(); iter++) {
    ENSURE(iter->GetType() == I3Particle::EMinus);
    ENSURE(iter->GetLocationType() == I3Particle::InActiveVolume);
    ENSURE(iter->GetMajorID() == p.GetMajorID());
    ENSURE(iter->GetMinorID() != p.GetMinorID() );

    // direction test
    ENSURE_DISTANCE(iter->GetZenith(), p.GetZenith(), 1e-7, "Direction of sub-cascade is not the same as original cascade");
    ENSURE_DISTANCE(iter->GetAzimuth(), p.GetAzimuth(), 1e-7, "Direction of sub-cascade is not the same as original cascade");
    totalEnergy += iter->GetEnergy();
  }
  // check energy
  ENSURE_DISTANCE(totalEnergy, p.GetEnergy(), 1e-7, "Total energy of sub-cascades is not the same as original cascade");
}

// creates particle of given energy and I3CascadeSplit instance
// passes particle to split algorithm, which uses parametrization or
// full simulation depending on energy, to split cascade into sub-cascades
// call Ensure which tests sub-cascade properties
void TestSplitting(double energy) {

  // cascade particle to split
  I3Particle p;
  p.SetType(I3Particle::EMinus);
  p.SetShape(I3Particle::Cascade);
  p.SetLocationType(I3Particle::InActiveVolume);
  // 23 deg zenith angle, 42 deg azimuth
  p.SetDir(23. * I3Units::deg, 42. * I3Units::deg);
  // energy 100 TeV 
  p.SetEnergy(energy);
  p.SetPos(0.,0.,0.);
  p.SetTime(0.);

  // random generator for I3CascadeSplit
  I3SPRNGRandomServicePtr random((new I3SPRNGRandomService(1,2,1)));
  // splitter instance
  I3CascadeSplit splitter(random);

  std::vector<I3Particle> subCascades = splitter.SplitCascade(p);
  // only run the following tests if a 'sane' energy
  // was passed.  we want to test pathological inputs
  // and ensure that I3CascadeSplit handles it well.
  // this means printing an error as opposed to aborting.
  if( std::isnormal( energy ) ){
    // make sure the splitting worked
    ENSURE(subCascades.size() > 0, "Cascade splitting failed");
    // run the tests on the sub-cascade vector
    Ensure(p, subCascades);
  }

}

/**********************************************
 * TEST DEFINITION SECTION
 **********************************************/

TEST_GROUP(I3CascadeSplit);

/* Test the cascade split algorithm
 *
 * Creates a cascade particle of 100 TeV
 * and passes it to the splitting algorithm
 *
 * The parameterization is used for splitting
 * Tests: 
 * starting point of sub-cascade is same as original
 * sum of sub-cascades energy is the same as original
 * direction of sub-cascades is the same as the original
 * 
 */
TEST(SplitUsingParametrization)
{
  TestSplitting(100 * I3Units::TeV);
}

/* Test the cascade split algorithm
 *
 * Creates a cascade particle with pathological
 * energies and passes it to the splitting algorithm
 * Tests: 
 * The splitter should handle this is in a dignified
 * manner, which means it should print an error and
 * tell the user that it's skipping over this particle.
 */
TEST(SplitUsingPathologicalInputs)
{
  TestSplitting(NAN);
  TestSplitting(std::numeric_limits<double>::infinity());
  TestSplitting(-std::numeric_limits<double>::infinity());
  TestSplitting(0.);
}

/* Test the cascade split algorithm
 *
 * Creates a cascade particle of 100 PeV
 * and passes it to the splitting algorithm
 *
 * The full simulation is used for splitting
 * Tests: 
 * starting point of sub-cascade is same as original
 * sum of sub-cascades energy is the same as original
 * direction of sub-cascades is the same as the original
 * 
 */
TEST(SplitUsingSimulation)
{
  TestSplitting(100 * I3Units::PeV);
}

TEST(MaxEnergy) {
  const double energy = 100 * I3Units::PeV;
  const double maxOutputEnergy = 100 * I3Units::TeV;
  
  // cascade particle to split
  I3Particle p;
  p.SetType(I3Particle::EMinus);
  p.SetShape(I3Particle::Cascade);
  p.SetLocationType(I3Particle::InActiveVolume);
  // 23 deg zenith angle, 42 deg azimuth
  p.SetDir(23. * I3Units::deg, 42. * I3Units::deg);
  // energy 100 TeV
  p.SetEnergy(energy);
  p.SetPos(0.,0.,0.);
  p.SetTime(0.);
  
  // random generator for I3CascadeSplit
  I3SPRNGRandomServicePtr random((new I3SPRNGRandomService(1,2,1)));
  // splitter instance
  I3CascadeSplit splitter(random);
  splitter.SetSegmentMaxEnergy(maxOutputEnergy);
  
  std::vector<I3Particle> subCascades = splitter.SplitCascade(p);
  // make sure the splitting worked
  ENSURE(subCascades.size() > 0, "Cascade splitting failed");
  // run the tests on the sub-cascade vector
  Ensure(p, subCascades);
  // check that all subcascades are below the requested threshold
  for(std::vector<I3Particle>::const_iterator it=subCascades.begin(), end=subCascades.end(); it!=end; it++){
    ENSURE(it->GetEnergy() <= maxOutputEnergy, "Cascades should be split to below the maximum allowed energy");
  }
  
}