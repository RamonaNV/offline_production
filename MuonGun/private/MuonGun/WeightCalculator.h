/** $Id: WeightCalculator.h 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 *$Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#ifndef I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED
#define I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <icetray/I3PointerTypedefs.h>

// only for GetMuonsAtSurface
#include <icetray/I3Frame.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>

#include <tableio/I3Converter.h>

namespace I3MuonGun {

I3_FORWARD_DECLARATION(Generator);
I3_FORWARD_DECLARATION(SamplingSurface);
I3_FORWARD_DECLARATION(Flux);
I3_FORWARD_DECLARATION(RadialDistribution);
I3_FORWARD_DECLARATION(EnergyDistribution);

std::vector<I3Particle>
GetMuonsAtSurface(I3FramePtr frame, I3Surfaces::SurfaceConstPtr surface);

/**
 * @brief Convenience struct to hold the components of a muon bundle flux model
 */
struct BundleModel {
	BundleModel(FluxPtr f, RadialDistributionPtr r, EnergyDistributionPtr e)
	    : flux(f), radius(r), energy(e) {}
	FluxPtr flux;
	RadialDistributionPtr radius;
	EnergyDistributionPtr energy;
};

/**
 * @brief Utility class to calculate weights for muon bundles
 */
class WeightCalculator {
public:
	/**
	 * @param[in] m Target model of the muon flux as a function of
	 *              angle, depth, multiplicity, radius, and energy
	 * @param[in] g Generation scheme according by which the
	 *              events were generated.
	 */
	WeightCalculator(const BundleModel &m, GenerationProbabilityPtr g)
	    : surface_(g->GetInjectionSurface()), flux_(m.flux), radius_(m.radius), energy_(m.energy), generator_(g) {}
	
	/**
	 * Calculate a weight corresponding to the frequency with which
	 * the given bundle configuration appears in the flux model.
	 *
	 * @param[in] axis   The bundle axis
	 * @param[in] bundle The radial offset and energy of each muon
	 *                   in the bundle
	 * @returns a weight in units of @f$ [s^{-1}] @f$
	 */
	double GetWeight(const I3Particle &axis, const BundleConfiguration &bundle) const;
	 
	SamplingSurfaceConstPtr GetSurface() { return surface_; }
	void SetSurface(SamplingSurfacePtr s) { surface_ = s; }

protected:
	WeightCalculator() {}
	SamplingSurfaceConstPtr surface_;
	FluxConstPtr flux_;
	RadialDistributionConstPtr radius_;
	EnergyDistributionConstPtr energy_;
	GenerationProbabilityConstPtr generator_;
};

class MuonBundleConverter : public I3ConverterImplementation<I3MCTree> {
public:
	MuonBundleConverter(size_t maxMultiplicity=25, SamplingSurfaceConstPtr surface=SamplingSurfaceConstPtr());
	I3TableRowDescriptionPtr CreateDescription(const I3MCTree&);
	size_t FillRows(const I3MCTree&, I3TableRowPtr);
private:
	SET_LOGGER("I3MuonGun::WeightCalculator");
	size_t maxMultiplicity_;
	SamplingSurfaceConstPtr surface_;
};

}

#endif // I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED
