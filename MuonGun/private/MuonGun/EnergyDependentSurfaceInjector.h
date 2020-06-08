/** $Id: EnergyDependentSurfaceInjector.h 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#ifndef I3MUONGUN_ENERGYDEPENDENTSURFACEINJECTOR_H_INCLUDED
#define I3MUONGUN_ENERGYDEPENDENTSURFACEINJECTOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <MuonGun/StaticSurfaceInjector.h>
#include <boost/function.hpp>

namespace I3MuonGun {

I3_FORWARD_DECLARATION(Cylinder);
I3_FORWARD_DECLARATION(SamplingSurface);
I3_FORWARD_DECLARATION(Flux);
I3_FORWARD_DECLARATION(OffsetPowerLaw);
I3_FORWARD_DECLARATION(RadialDistribution);

class SurfaceScalingFunction {
public:
	virtual ~SurfaceScalingFunction();
	
	/** @brief Propose a target surface for the given energy */
	virtual SamplingSurfacePtr GetSurface(double energy) const = 0;
	
	/** @brief Compare for equality */
	virtual bool operator==(const SurfaceScalingFunction&) const = 0;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(SurfaceScalingFunction);

class ConstantSurfaceScalingFunction : public SurfaceScalingFunction {
public:
	ConstantSurfaceScalingFunction(SamplingSurfacePtr surface);
	virtual ~ConstantSurfaceScalingFunction();
	
	virtual SamplingSurfacePtr GetSurface(double energy) const override;
	virtual bool operator==(const SurfaceScalingFunction&) const override;
private:
	ConstantSurfaceScalingFunction();
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	SamplingSurfacePtr surface_;
};

I3_POINTER_TYPEDEFS(ConstantSurfaceScalingFunction);

class BasicSurfaceScalingFunction : public SurfaceScalingFunction {
public:
	BasicSurfaceScalingFunction();
	virtual ~BasicSurfaceScalingFunction();
	
	virtual SamplingSurfacePtr GetSurface(double energy) const override;
	virtual bool operator==(const SurfaceScalingFunction&) const override;
	
	void SetCapScaling(double energyScale, double scale, double offset, double power);
	void SetSideScaling(double energyScale, double scale, double offset, double power);
	void SetRadiusBounds(double rmin, double rmax);
	void SetZBounds(double zmin, double zmax);
	
	void SetCenterMin(double x, double y);
	void SetCenterMax(double x, double y);
	
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	double GetMargin(double logenergy, double scale, double offset, double power) const;
	
	typedef std::pair<double, double> pair;
	pair scale_, energyScale_, offset_, power_;
	pair rBounds_;
	pair zBounds_;
	std::pair<pair, pair> centerBounds_;
};

I3_POINTER_TYPEDEFS(BasicSurfaceScalingFunction);

/**
 * @brief A rejection-sampling Generator with energy-dependent sampling surface
 *
 * EnergyDependentSurfaceInjector samples bundle impact points, angles, multiplicities,
 * and radial distributions at their natural frequencies, but scales the sampling
 * surface based on the highest-energy muon in the bundle: dim, low-energy muons are
 * targeted only at a small inner surface, while the surface scales up to full size
 * for potentially bright muons. This technique can be used to efficiently simulate
 * background for an event selection that requires a thick veto for dim events
 * (where the rates are also highest) but becomes more accepting for bright events.
 */
class EnergyDependentSurfaceInjector : public StaticSurfaceInjector {
public:
	EnergyDependentSurfaceInjector(CylinderPtr surface=CylinderPtr(), FluxPtr flux=FluxPtr(), boost::shared_ptr<OffsetPowerLaw> energies=boost::shared_ptr<OffsetPowerLaw>(),
	    RadialDistributionPtr radius=RadialDistributionPtr(),
	    SurfaceScalingFunctionPtr scaling=boost::make_shared<BasicSurfaceScalingFunction>());

	// GenerationProbability interface
	virtual double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const override;
	virtual GenerationProbabilityPtr Clone() const override;
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const override;
	
	// Generator interface
	void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	
	SurfaceScalingFunctionPtr GetScaling() const { return scalingFunction_; }
	void SetScaling(SurfaceScalingFunctionPtr &f) { scalingFunction_ = f; }
	
	/** 
	 * Scale the sampling cylinder to a size appropriate for
	 * the given maximum muon energy
	 * by the given energy. This is not necessarily the fastest.
	 */
	SamplingSurfacePtr GetTargetSurface(double energy) const;
	/** 
	 * Integrate the flux to get the total rate on the surface.
	 * This is not necessarily the fastest.
	 */
	double GetTotalRate(SamplingSurfaceConstPtr surface) const;
private:
	
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	SurfaceScalingFunctionPtr scalingFunction_;
};

}

I3_CLASS_VERSION(I3MuonGun::SurfaceScalingFunction, 0);
I3_CLASS_VERSION(I3MuonGun::ConstantSurfaceScalingFunction, 0);
I3_CLASS_VERSION(I3MuonGun::BasicSurfaceScalingFunction, 0);
I3_CLASS_VERSION(I3MuonGun::EnergyDependentSurfaceInjector, 0);

#endif
