/** $Id: RadialDistribution.h 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#ifndef I3MUONGUN_RADIALDISTRIBUTION_H
#define I3MUONGUN_RADIALDISTRIBUTION_H

#include <MuonGun/SplineTable.h>
#include <icetray/I3PointerTypedefs.h>

class I3Position;
class I3RandomService;

namespace I3MuonGun {

/**
 * @brief The distribution of distance between a muon and the bundle axis
 */
class RadialDistribution {
public:
	virtual ~RadialDistribution();
	/**
	 * @brief Calculate the probability of obtaining the given radial offset
	 *
	 * @param[in] depth        vertical depth in km
	 * @param[in] cos_theta    cosine of zenith angle
	 * @param[in] multiplicity number of muons in the bundle
	 * @param[in] radius       distance to bundle axis
	 * @returns a properly normalized propability @f$ dP/dr \,\, [m^{-1}] @f$
	 */
	double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	
	/**
	 * @brief Calculate the logarithm of the probability of obtaining
	 * the given radial offset
	 *
	 * @see operator()
	 */
	virtual double GetLog(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const = 0;
	
	/**
	 * @brief Draw a sample from the distribution of radii
	 *
	 * @param[in] rng          a random number generator
	 * @param[in] depth        vertical depth in km
	 * @param[in] cos_theta    cosine of zenith angle
	 * @param[in] multiplicity number of muons in the bundle
	 * @returns a radius in m
	 */
	virtual double Generate(I3RandomService &rng, double depth, double cos_theta,
	    unsigned multiplicity) const = 0;
	    
	
	virtual bool operator==(const RadialDistribution&) const = 0;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(RadialDistribution);

/**
 * @brief Radial distribution according to Becherini et al.
 */
class BMSSRadialDistribution : public RadialDistribution {
public:
	BMSSRadialDistribution();
	double GetLog(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	double Generate(I3RandomService &rng, double depth, double cos_theta,
	    unsigned multiplicity) const;
	
	virtual bool operator==(const RadialDistribution&) const override;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	double GetMeanRadius(double, double, unsigned) const;
	double GetShapeParameter(double, double, unsigned) const;
	double GetGenerationProbability(double, double, double) const;
	double rho0a_, rho0b_, rho1_, theta0_, f_, alpha0a_, alpha0b_, alpha1a_, alpha1b_;
	double rmax_;
};

I3_POINTER_TYPEDEFS(BMSSRadialDistribution);

/**
 * @brief Radial distribution fit to a tensor-product B-spline surface
 *
 * The surface is fit to @f$ d \log{P} / d{r^2} @f$ to remove the factor
 * of differential area implicit in @f$ dP / dr @f$
 */
class SplineRadialDistribution : public RadialDistribution {
public:
	SplineRadialDistribution(const std::string&);
	double GetLog(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	double Generate(I3RandomService &rng, double depth, double cos_theta,
	    unsigned multiplicity) const;
	
	virtual bool operator==(const RadialDistribution&) const override;
private:
	SplineRadialDistribution() {}
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	SplineTable spline_;
};

}

I3_CLASS_VERSION(I3MuonGun::RadialDistribution, 0);
I3_CLASS_VERSION(I3MuonGun::SplineRadialDistribution, 0);

#endif // I3MUONGUN_RADIALDISTRIBUTION_H
