/** $Id: NaturalRateInjector.h 159352 2017-11-07 22:43:41Z juancarlos $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 159352 $
 * $Date: 2017-11-07 15:43:41 -0700 (Tue, 07 Nov 2017) $
 */

#ifndef MUONGUN_NATURALRATEINJECTOR_H_INCLUDED
#define MUONGUN_NATURALRATEINJECTOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <MuonGun/SamplingSurface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <dataclasses/physics/I3Particle.h>

#include <boost/tuple/tuple.hpp>

namespace I3MuonGun {

/**
 * @brief A natural-rate Generator
 *
 * NaturalRateInjector samples bundle impact points, angles, multiplicities,
 * and radius/energy distributions at their natural frequencies on a fixed
 * surface using a brain-dead acceptance/rejection technique. This is
 * to MUPAGE, except that the flux parameterizations can be swapped out and the
 * events can be combined with biased generation and reweighted to a different
 * parameterization after the fact. 
 */
class NaturalRateInjector : public Generator {
public:
	NaturalRateInjector();
	NaturalRateInjector(SamplingSurfacePtr surface, FluxPtr flux,
	    EnergyDistributionPtr edist);
	
	// Generator Interface
	virtual void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	virtual GenerationProbabilityPtr Clone() const;
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const;
	virtual double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
	virtual SamplingSurfaceConstPtr GetInjectionSurface() const { return surface_; }
	
	void SetSurface(SamplingSurfacePtr p);
	SamplingSurfacePtr GetSurface() { return surface_; }
	
	void SetFlux(FluxPtr p);
	FluxPtr GetFlux() { return flux_; }
	
	void SetEnergyDistribution(EnergyDistributionPtr e);
	EnergyDistributionPtr GetEnergyDistribution() { return energyDistribution_; }
	
	/**
	 * Integrate the configured flux over the sampling surface, summing over
	 * all allowed multiplicities.
	 *
	 * @returns a rate in units of @f$ [s^{-1}] @f$
	 */
	double GetTotalRate() const;

protected:
	/**
	 * Draw a sample from the distribution of shower impact points
	 *
	 * The shower axis and multiplcity are filled into the *axis*
	 */
	void GenerateAxis(I3RandomService &rng, std::pair<I3Particle, unsigned> &axis) const;
	/**
	 * Distribute the given number of muons in the transverse plane
	 * and draw an energy for each
	 */
	void FillMCTree(I3RandomService &rng, const std::pair<I3Particle, unsigned> &axis, I3MCTree &, BundleConfiguration &) const;

private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);

protected:
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	EnergyDistributionPtr energyDistribution_;
	
	mutable double totalRate_;

};

}

I3_CLASS_VERSION(I3MuonGun::NaturalRateInjector, 0);

#endif

