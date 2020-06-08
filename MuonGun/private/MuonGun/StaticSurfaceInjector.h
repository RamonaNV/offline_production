/** $Id: StaticSurfaceInjector.h 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#ifndef MUONGUN_CANCAN_H_INCLUDED
#define MUONGUN_CANCAN_H_INCLUDED

#include <MuonGun/Generator.h>
#include <MuonGun/SamplingSurface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <dataclasses/physics/I3Particle.h>

#include <boost/tuple/tuple.hpp>

namespace I3MuonGun {

/**
 * @brief A simple rejection-sampling Generator
 *
 * StaticSurfaceInjector samples bundle impact points, angles, multiplicities,
 * and radial distributions at their natural frequencies on a fixed surface
 * using a brain-dead acceptance/rejection technique. Energies, on the other hand,
 * are sampled from an OffsetPowerLaw, which both boosts efficiency and allows
 * a measure of control over the generated energy distribution.
 */
class StaticSurfaceInjector : public Generator {
public:
	StaticSurfaceInjector();
	StaticSurfaceInjector(SamplingSurfacePtr surface, FluxPtr flux,
	    boost::shared_ptr<OffsetPowerLaw> edist, RadialDistributionPtr rdist);
	
	// Generator Interface
	virtual void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const override;
	virtual GenerationProbabilityPtr Clone() const override;
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const override;
	virtual double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const override;
	virtual SamplingSurfaceConstPtr GetInjectionSurface() const override { return surface_; }
	
	void SetSurface(SamplingSurfacePtr p);
	SamplingSurfacePtr GetSurface() { return surface_; }
	
	void SetFlux(FluxPtr p);
	FluxPtr GetFlux() { return flux_; }
	
	void SetRadialDistribution(RadialDistributionPtr r) { radialDistribution_ = r; }
	RadialDistributionPtr GetRadialDistribution() { return radialDistribution_; }
	
	void SetEnergyDistribution(boost::shared_ptr<OffsetPowerLaw> e) { energyGenerator_ = e; }
	boost::shared_ptr<OffsetPowerLaw> GetEnergyDistribution() { return energyGenerator_; }
	
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
	
	void CalculateMaxFlux();
	
	/**
	 * Get the normalization term for relative weighting of zenith
	 * angles and multiplicities by integrating the flux at the top
	 * of the cylinder over all zenith angles and summing over all
	 * allowed multiplicities.
	 */
	double GetZenithNorm() const;

private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);

protected:
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	RadialDistributionPtr radialDistribution_;
	
	double maxFlux_;
	mutable double totalRate_, zenithNorm_;

};

}

I3_CLASS_VERSION(I3MuonGun::StaticSurfaceInjector, 1);

#endif

