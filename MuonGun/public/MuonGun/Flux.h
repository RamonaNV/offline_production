/** $Id: Flux.h 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#ifndef I3MUONGUN_FLUX_H_INCLUDED
#define I3MUONGUN_FLUX_H_INCLUDED

#include <MuonGun/SplineTable.h>
#include <icetray/I3PointerTypedefs.h>

namespace I3MuonGun {

/**
 * @brief A parameterization of the total muon-bundle flux
 * 
 * The total flux (in units of @f$ [ m^{-2} sr^{-2} s^{-1} ]@f$)
 * is parameterized in terms of vertical depth @f$ km @f$, the
 * cosine of the zenith angle, and bundle multiplicity.
 */
class Flux {
public:
	Flux();
	virtual ~Flux();
	typedef double result_type;
	double operator()(double depth, double cos_theta, unsigned multiplicity) const;
	virtual double GetLog(double depth, double cos_theta, unsigned multiplicity) const = 0;

	unsigned GetMaxMultiplicity() const { return maxMultiplicity_; }
	unsigned GetMinMultiplicity() const { return minMultiplicity_; }
	
	void SetMaxMultiplicity(unsigned m) { maxMultiplicity_ = m; }
	void SetMinMultiplicity(unsigned m) { minMultiplicity_ = m; }
	
	virtual bool operator==(const Flux&) const;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	unsigned minMultiplicity_, maxMultiplicity_;
};

I3_POINTER_TYPEDEFS(Flux);

/**
 * @brief Total flux according to Becherini et al.
 */
class BMSSFlux : public Flux {
public:
	BMSSFlux();
	double GetLog(double depth, double cos_theta, unsigned multiplicity) const;
	
	virtual bool operator==(const Flux&) const override;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	double k0a_, k0b_, k1a_, k1b_;
	double v0a_, v0b_, v0c_, v1a_, v1b_;
};

/**
 * @brief Total flux, fit to a tensor-product B-spline surface
 */
class SplineFlux : public Flux {
public:
	SplineFlux(const std::string &singles, const std::string &bundles);
	double GetLog(double depth, double cos_theta, unsigned multiplicity) const;
	
	virtual bool operator==(const Flux&) const override;
private:
	SplineFlux() {}
	
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	SplineTable singles_;
	SplineTable bundles_;
};

}

#endif // I3MUONGUN_FLUX_H_INCLUDED
