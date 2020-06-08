#include <I3Test.h>

#include <vector>
#include <limits>
#include "phys-services/I3SPRNGRandomService.h"
#include "icetray/I3Units.h"
#include "dataclasses/physics/I3Particle.h"
#include "cmc/I3CascadeMCService.h"
#include "cmc/I3MetropolisHastings.h"
#include "cmc/I3CascadeSimulation.h"
#include "cmc/I3CascadeSimulationCrossSections.h"

TEST_GROUP(Repeatability);

TEST(MetropolisSampler)
{
	I3RandomServicePtr rng(new I3SPRNGRandomService(0,1,0));
	
	beta_dist::alpha = 0.5;
	beta_dist::beta = 0.5;
	I3MetropolisHastings sampler(rng,
	    I3CascadeSimulationCrossSections::BremsstrahlungCrossSection,
	    &beta_dist::rvs, &beta_dist::pdf, "testy");
	
	std::vector<double> samples;
	
	I3FrameObjectPtr state = rng->GetState();
	sampler.BurnIn(I3Units::TeV, .1, 100);
	for (int i=0; i < 100; i++)
		samples.push_back(sampler.Sample(std::pow(10, rng->Uniform(3, 7))));
	
	for (int j=0; j < 10; j++) {
		rng->RestoreState(state);
		sampler.BurnIn(I3Units::TeV, .1, 100);
		for (int i=0; i < 100; i++) {
			double sample = sampler.Sample(std::pow(10, rng->Uniform(3, 7)));
			ENSURE_EQUAL(sample, samples[i]);
		}
	}
}

TEST(CascadeSimulation)
{
	I3RandomServicePtr rng(new I3SPRNGRandomService(0,1,0));
	
	I3CascadeSimulation sim(rng);
	// Activate per-call burn-in
	sim.SetRandomNumberGenerator(rng);
	
	double energy = 1e7;
	I3FrameObjectPtr state = rng->GetState();
	sim.Simulate(I3Particle::EMinus, energy);
	std::vector<double> eloss = sim.GetEnergyLossProfile();
	
	for (int j=0; j < 10; j++) {
		rng->RestoreState(state);
		sim.Simulate(I3Particle::EMinus, energy);
		std::vector<double> eloss_comp = sim.GetEnergyLossProfile();
		ENSURE_EQUAL(eloss_comp.size(), eloss.size());
		for (unsigned i=0; i < eloss.size(); i++)
			ENSURE_EQUAL(eloss_comp[i], eloss[i]);
	}
}

TEST(CascadeMCService)
{
	I3RandomServicePtr rng(new I3SPRNGRandomService(2,10000,1));
	
	I3CascadeMCService sim(rng);
	// Activate per-call burn-in
	sim.SetRandomNumberGenerator(rng);
	
	I3Particle p;
	p.SetPos(0,0,0);
	p.SetDir(0,0);
	p.SetTime(0);
	p.SetType(I3Particle::DeltaE);
	p.SetLocationType(I3Particle::InIce);
	p.SetEnergy(1e7);
	
	std::vector<I3Particle> daughters;
	I3FrameObjectPtr state = rng->GetState();
	sim.Simulate(p, daughters);
	ENSURE(daughters.size() > 0);
	
	for (int j=0; j < 10; j++) {
		std::vector<I3Particle> daughters_comp;
		rng->RestoreState(state);
		sim.Simulate(p, daughters_comp);
		ENSURE_EQUAL(daughters_comp.size(), daughters.size());
		for (unsigned i=0; i < daughters.size(); i++) {
			ENSURE_EQUAL(daughters_comp[i].GetTime(), daughters[i].GetTime());
			ENSURE_EQUAL(daughters_comp[i].GetEnergy(), daughters[i].GetEnergy());
		}
	}
}
