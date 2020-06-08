/** $Id: SamplingSurface.h 146654 2016-06-01 18:31:33Z cweaver $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 146654 $
 * $Date: 2016-06-01 12:31:33 -0600 (Wed, 01 Jun 2016) $
 */

#ifndef I3MUONGUN_SAMPLINGSURFACE_H_INCLUDED
#define I3MUONGUN_SAMPLINGSURFACE_H_INCLUDED

#include <phys-services/surfaces/SamplingSurface.h>
#include <boost/function.hpp>

class I3Direction;
class I3RandomService;

namespace I3MuonGun {

/**
 * @brief A surface upon which muon bundles can be generated
 *
 * SamplingSurface knows how to calculate its projected area and integrate
 * a flux over its surface. It is assumed to be azimuthally symmetric, but
 * its projected area may vary with zenith angle.
 */
class SamplingSurface : public I3Surfaces::SamplingSurface {
public:
	virtual ~SamplingSurface();
	/** Get the integral of area*solid_angle over the given cos(theta) range */
	virtual double GetAcceptance(double cosMin=0, double cosMax=1) const = 0;
	/** Get the minimum vertical depth the surface occupies */
	virtual double GetMinDepth() const = 0;

	/** 
	 * Integrate a flux (defined in terms of dN/dOmega(depth [km], cos(theta)))
	 * over the outer surface
	 */
	virtual double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const = 0;

	virtual bool operator==(const SamplingSurface&) const = 0;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(SamplingSurface);

}

I3_CLASS_VERSION(I3MuonGun::SamplingSurface, 0);

#endif // I3MUONGUN_SAMPLINGSURFACE_H_INCLUDED
