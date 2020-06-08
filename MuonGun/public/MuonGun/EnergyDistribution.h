/** $Id: EnergyDistribution.h 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#ifndef I3MUONGUN_ENERGYDISTRIBUTION_H_INCLUDED
#define I3MUONGUN_ENERGYDISTRIBUTION_H_INCLUDED

#include <icetray/I3Units.h>
#include <icetray/serialization.h>
#include <icetray/I3PointerTypedefs.h>
#include <MuonGun/SplineTable.h>

class I3RandomService;

namespace I3MuonGun {

class RadialDistribution;
class OffsetPowerLaw;

/**
 * @brief Normalized distribution of energies within a bundle
 *
 * The energy distribution @f$ dP/dE \,\, [GeV^{-1}] @f$ is
 * parameterized in terms of vertical depth @f$ [km] @f$,
 * cosine of the zenith angle, bundle multiplicity, perpendicular
 * distance from the bundle axis @f$ [m] @f$, and energy  @f$ [GeV] @f$.
 *
 * Single muons are assumed to be colinear with the bundle axis.
 */
class EnergyDistribution {
public:
	EnergyDistribution() : minLog_(std::log(I3Units::GeV)), maxLog_(std::log(I3Units::PeV)) {}
	virtual ~EnergyDistribution();
	typedef double result_type;
	double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius, double energy) const;
	double Integrate(double depth, double cos_theta, unsigned multiplicity,
	    double r_min, double r_max, double e_min, double e_max) const;

	/// A tag type to enforce units
	struct log_value {
		explicit log_value(double v) : value_(v) {}
		operator double() const { return value_; }
		double value_;
	};
	virtual double GetLog(double depth, double cos_theta,
	    unsigned multiplicity, double radius, log_value log_energy) const = 0;

	/// Sample *samples* (radius, energy) pairs
	virtual std::vector<std::pair<double,double> > Generate(I3RandomService &rng,
	    double depth, double cos_theta, unsigned multiplicity, unsigned samples) const = 0;
	
	double GetMax() const { return std::exp(maxLog_); }
	double GetMin() const { return std::exp(minLog_); }
	void SetMax(double v) { maxLog_ = std::log(v); }
	void SetMin(double v) { minLog_ = std::log(v); }
	
	virtual double GetMaxRadius() const = 0;
	
	virtual bool operator==(const EnergyDistribution&) const = 0;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
protected:
	double minLog_, maxLog_;
};

I3_POINTER_TYPEDEFS(EnergyDistribution);

/**
 * @brief Energy distribution fit to a tensor-product B-spline surface
 *
 * The surface is fit to @f$ d \log{P} / d\log(E) @f$, which is nearly
 * polynomial and thus easier to represent with low-order splines.
 */
class SplineEnergyDistribution : public EnergyDistribution {
public:
	SplineEnergyDistribution(const std::string &singles, const std::string &bundles);
	double GetLog(double depth, double cos_theta, 
	    unsigned multiplicity, double radius, log_value log_energy) const;
	std::vector<std::pair<double,double> > Generate(I3RandomService &rng,
	    double depth, double cos_theta, unsigned multiplicity, unsigned samples) const;
	virtual double GetMaxRadius() const override;
	virtual bool operator==(const EnergyDistribution&) const override;
private:
	SplineEnergyDistribution() {}
	
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);

	SplineTable singles_;
	SplineTable bundles_;
};

class BMSSEnergyDistribution : public EnergyDistribution {
public:
	BMSSEnergyDistribution();
	double GetLog(double depth, double cos_theta, 
	    unsigned multiplicity, double radius, log_value log_energy) const;
	std::vector<std::pair<double,double> > Generate(I3RandomService &rng,
	    double depth, double cos_theta, unsigned multiplicity, unsigned samples) const;
	
	virtual double GetMaxRadius() const override;
	
	virtual bool operator==(const EnergyDistribution&) const override;
	
	OffsetPowerLaw GetSpectrum(double depth, double cos_theta, unsigned m, double r) const;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	// Single-muon energy distribution
	double beta_, g0_, g1_, e0a_, e0b_, e1a_, e1b_;
	// Bundle energy distribution
	double a0_, a1_, b0a_, b0b_, b1a_, b1b_, q0_, q1_,
	    c0a_, c0b_, c1_, d0a_, d0b_, d1a_, d1b_;
};

/**
 * @brief An approximate form for the underground muon energy spectrum
 *
 * The deep underground muon spectrum very nearly follows a power law with
 * a break at low energies caused by pile-up of minimum-ionizing muons.
 * While the details of the spectrum depend on depth, zenith angle, distance
 * from the bundle axis, etc., this approximate form can still be used to
 * efficiently generate plausible energies that can later be weighted to
 * reflect their frequency under a more precise model.
 */
class OffsetPowerLaw {
public:
	OffsetPowerLaw();
	/**
	 * Create an offset power-law energy distribution of the form
	 * @f$ dP/dE_{\mu} \propto (E_{\mu} + b)^{-\gamma} @f$ 
	 *
	 * @param[in] gamma  The power law index (positive)
	 * @param[in] offset The offset in @f$ (E_{\mu} + b) @f$
	 * @param[in] emin   minimum generated energy
	 * @param[in] emax   maximum generated energy
	 */
	OffsetPowerLaw(double gamma, double offset, double emin, double emax);
	typedef double result_type;
	/** Calculate the probability that the given energy was generated */
	double operator()(double energy) const;
	double GetLog(double energy) const;
	double GetLog(EnergyDistribution::log_value log_energy) const;
	/** Draw an energy from the distribution */
	double Generate(I3RandomService &rng) const;
	double InverseSurvivalFunction(double p) const;
	
	const double GetMin() const { return emin_; }
	const double GetMax() const { return emax_; }
	
	bool operator==(const OffsetPowerLaw &other) const;
private:
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	double gamma_, offset_;
	double emin_, emax_;
	double nmin_, nmax_, norm_, lognorm_;
};

}

I3_CLASS_VERSION(I3MuonGun::EnergyDistribution, 0);
I3_CLASS_VERSION(I3MuonGun::SplineEnergyDistribution, 0);
I3_CLASS_VERSION(I3MuonGun::BMSSEnergyDistribution, 0);
I3_CLASS_VERSION(I3MuonGun::OffsetPowerLaw, 0);

#endif
