#include <I3Test.h>

#include <vector>
#include <cmath>  
#include "phys-services/I3SPRNGRandomService.h"
#include "phys-services/I3SPRNGRandomServiceFactory.h"
#include "dataclasses/physics/I3Particle.h"
#include "icetray/I3Units.h" 
#include "cmc/I3CascadeMuonSplit.h"


// checks the direction, position and total energy of sub-cascade list
void MuonEnsure(I3Particle p, double energy, 
		std::vector<I3Particle> subCascades) {

  std::vector<I3Particle>::const_iterator iter = subCascades.begin();


  // check direction and sum energies

  double totalEnergy = 0.0;
  for(; iter != subCascades.end(); iter++) {

    ENSURE(iter->GetType() == I3Particle::MuPlus);
    ENSURE(iter->GetLocationType() == I3Particle::InIce);
    ENSURE(iter->GetMajorID() == p.GetMajorID());
    ENSURE(iter->GetMinorID() != p.GetMinorID() );
    ENSURE(std::isnan(iter->GetLength() )); 
    // check location, this is a float comparison
    ENSURE_DISTANCE(iter->GetX(), p.GetX(), 1e-7, 
		    "Start location of sub-cascades is "	\
		    "not the same as the original cascade");
    ENSURE_DISTANCE(iter->GetY(), p.GetY(), 1e-7, 
		    "Start location of sub-cascades is "	\
		    "not the same as the original cascade");
    ENSURE_DISTANCE(iter->GetZ(), p.GetZ(), 1e-7, 
		    "Start location of sub-cascades is "	\
		    "not the same as the original cascade");

    // direction test
    ENSURE_DISTANCE(iter->GetZenith(), p.GetZenith(), 1e-7, 
		    "Direction of sub-cascade is not the "\
		    "same as original cascade");
    ENSURE_DISTANCE(iter->GetAzimuth(), p.GetAzimuth(), 1e-7, 
		    "Direction of sub-cascade is not the "\
		    "same as original cascade");
    totalEnergy += iter->GetEnergy();
  }
  // check energy
  ENSURE_DISTANCE(totalEnergy + p.GetEnergy(), energy, 1e-7, 
		  "Total energy muons plus rest cascade is "\
		  "not the same as original cascade");
}

// creates particle of given energy and I3CascadeMuonSplit instance
// passes particle to split algorithm, to generate muons and reduce 
// the cascade energy. 
// call Ensure which tests sub-cascade properties
void TestMuonSplitting(double energy) {

  // cascade particle to split
  I3Particle p;
  // random generator for I3CascadeSplit
  I3SPRNGRandomServicePtr random((new I3SPRNGRandomService(1,2,1)));
  // splitter instance
  I3CascadeMuonSplit splitter = I3CascadeMuonSplit(random);
  std::vector<I3Particle> subCascades;

  p.SetType(I3Particle::PPlus);
  p.SetShape(I3Particle::Cascade);
  // 23 deg zenith angle, 42 deg azimuth
  p.SetDir(23. * I3Units::deg, 42. * I3Units::deg);
  // energy 100 TeV 
  p.SetEnergy(energy);
  p.SetPos(0.,0.,0.);
  p.SetTime(0.);
  subCascades = splitter.GenerateMuons(p);

  // run the following tests if the energy is "normal"
  // otherwise cmc should pass over this particle gracefully
  // while printing an error message.
  if( std::isnormal( energy ) ){
    // make sure the splitting worked
    ENSURE(subCascades.size() > 0, "Cascade muon splitting failed");
    // run the tests on the sub-cascade vector
    MuonEnsure(p, energy, subCascades);
  }
}

/**********************************************
 * TEST DEFINITION SECTION
 **********************************************/

TEST_GROUP(I3CascadeMuonSplit);

/* Test the cascade muon split algorithm
 *
 * Creates a (hadronic) cascade particle of 100 TeV
 * and passes it to the muon splitting algorithm
 *
 * Tests: 
 * starting point of muons is same as original
 * sum of muon energy plus rest cascade is the same as original
 * direction of muons is the same as the original
 * 
 */
TEST(Split)
{
  TestMuonSplitting(100 * I3Units::TeV);
}

TEST(SplitUsingPathologicalInputs)
{
  TestMuonSplitting(NAN);
  TestMuonSplitting(std::numeric_limits<double>::infinity());
  TestMuonSplitting(-std::numeric_limits<double>::infinity());
  TestMuonSplitting(0.);
}



