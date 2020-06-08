/**
 * brief Implementation I3CascadeSplit class
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: I3CascadeSplit.cxx 150508 2016-10-04 14:44:25Z jvansanten $
 *
 * @version $Revision: 150508 $
 * @date $LastChangedDate: 2016-10-04 08:44:25 -0600 (Tue, 04 Oct 2016) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: jvansanten $
 */

//Standard C/C++ includes
#include <cmath>
#include <vector>
#include <algorithm>

// local includes
#include "I3CascadeDevelopment.h"
#include "I3CascadeSplit.h"

//IceTray includes
#include "dataclasses/I3Constants.h"
#include "icetray/I3Units.h"
//#include "phys-services/I3TRandomService.h"

using namespace std;
// Class constants
const double I3CascadeSplit::DEFAULT_SIMULATION_THRESHOLD = 1 * I3Units::PeV;
const int I3CascadeSplit::DEFAULT_STEP_WIDTH = 3;  // units are radiation length
const double I3CascadeSplit::DEFAULT_MAX_ENERGY = 1 * I3Units::PeV;

// constructor
I3CascadeSplit::I3CascadeSplit(I3RandomServicePtr random) :
  cascadeSimulation_(random),
  enableSimulation_(true),
  simulationThreshold_(I3CascadeSplit::DEFAULT_SIMULATION_THRESHOLD),
  stepWidth_(I3CascadeSplit::DEFAULT_STEP_WIDTH),
  segmentMaxEnergy_(DEFAULT_MAX_ENERGY)
{
	SetRandomNumberGenerator(random);
}

I3CascadeSplit::~I3CascadeSplit() {}

void
I3CascadeSplit::SetRandomNumberGenerator(I3RandomServicePtr r)
{
	I3CascadeMCCommon::Sampler::SetRandomNumberGenerator(r);
	cascadeSimulation_.SetRandomNumberGenerator(r);
	cascadeParametrization_.SetRandomNumberGenerator(r);
}

// Implementation of the cascade split method.
vector<I3Particle> I3CascadeSplit::SplitCascade(I3Particle &cascade) {
  log_trace("Entering SplitCascade");
  log_debug("Cascade: type %s at (%.2f, %.2f, %.2f) with energy (%.3e)", 
            cascade.GetTypeString().c_str(), 
            cascade.GetX(), 
            cascade.GetY(),
            cascade.GetZ(), 
            cascade.GetEnergy());

  if( !std::isnormal(cascade.GetEnergy()) ){
    log_error("The cascade has a pathological energy = %e",
	      cascade.GetEnergy());
    log_error("Skipping this cascade.");
    return vector<I3Particle>();
  }

  // simulate the shower (calculation of parametrization or full simulation)
  // if cascade simulation is enabled and energy is larger than the threshold, 
  // simulate the shower otherwise use the parametrization
  I3CascadeDevelopment* cascadeDevelopment;

  if (enableSimulation_ && (cascade.GetEnergy() >= simulationThreshold_)) {
    log_debug("Simulate shower (energy: %.3e)", cascade.GetEnergy());
    cascadeDevelopment = &cascadeSimulation_;
  } else {
    log_debug("Parametrize shower (energy: %.3e)", cascade.GetEnergy());
    cascadeDevelopment = &cascadeParametrization_;
  }

  // the simulation is only one dimensional, no position and direction information needed
  // compute the longitudinal development of the shower
  cascadeDevelopment->Simulate(cascade.GetType(), cascade.GetEnergy());

  // get the length in units of radiation lenght
  double length = cascadeDevelopment->GetLength();
  log_debug("Shower length in radiation length %f", length);

  // init the vector for the sub cascades
  vector<I3Particle> shower;
  
  // loop until the end of the shower and create sub-cascades
  for(int i = 0; i < length; i = i+stepWidth_) {
    log_trace("Creating sub-cascade #%i", i);

    // get the energy loss between this step and the next
    double energy = cascadeDevelopment->GetEnergyDeposit(i, i+stepWidth_);
    log_trace("Fractional energy %f", energy/cascade.GetEnergy());

    while(energy>=0){
      // create a new particle i steps away from the original
      double energyChunk = std::min(energy, segmentMaxEnergy_);
      I3Particle subCascade = CreateParticle(cascade, energyChunk, i, stepWidth_);
      subCascade.SetShape(I3Particle::CascadeSegment);
      subCascade.SetLength(stepWidth_ * I3CascadeDevelopment::RADIATION_LENGTH * I3Units::m);
      shower.push_back(subCascade);
      log_debug("Added Subcascade X:%.2f Y:%.2f Z:%.2f T:%.2f E:%.3e",
                subCascade.GetX(), subCascade.GetY(), subCascade.GetZ(),
                subCascade.GetTime(), subCascade.GetEnergy());
      energy -= energyChunk;
      if(energy==0) // special case for zero energy cascades
        break;
    }
  }
  log_trace("Leaving SplitCascade");
  return shower;
}

// creates a EMinus particle with given energy, lenght and i*radiationlength distance to the orignal position
// in the direction of the original particle
I3Particle I3CascadeSplit::CreateParticle(I3Particle &particle, double energy, int step, int length) {
  
  I3Particle sub_particle;
  double i = step + 0.5;
  sub_particle.SetLocationType(particle.GetLocationType());
  sub_particle.SetDir(particle.GetDir());
  sub_particle.SetType(I3Particle::EMinus); // sub cascade is em shower
  sub_particle.SetEnergy(energy);
  sub_particle.SetTime(particle.GetTime() + i * I3CascadeDevelopment::RADIATION_LENGTH * I3Units::m / I3Constants::c);
  I3Position position = sub_particle.GetPos();
  position.SetX(particle.GetX() + particle.GetDir().GetX() * \
                i * I3CascadeDevelopment::RADIATION_LENGTH * I3Units::m);
  position.SetY(particle.GetY() + particle.GetDir().GetY() * \
                i * I3CascadeDevelopment::RADIATION_LENGTH * I3Units::m);
  position.SetZ(particle.GetZ() + particle.GetDir().GetZ() * \
                i * I3CascadeDevelopment::RADIATION_LENGTH * I3Units::m);
  sub_particle.SetPos(position);

  return sub_particle;
}


    
    

  


