/** $Id: EnergyDistribution.cxx 165001 2018-08-27 09:22:30Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 165001 $
 * $Date: 2018-08-27 03:22:30 -0600 (Mon, 27 Aug 2018) $
 */

#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/RadialDistribution.h>
#include <phys-services/I3RandomService.h>
#include <MuonGun/EnsembleSampler.h>

#include <boost/bind.hpp>

namespace I3MuonGun {

EnergyDistribution::~EnergyDistribution() {};

double
EnergyDistribution::operator()(double d, double ct, 
    unsigned m, double r, double e) const
{
	return std::exp(GetLog(d, ct, m, r, log_value(std::log(e))));
}

double
EnergyDistribution::Integrate(double d, double ct, 
    unsigned m, double r_min, double r_max, double e_min, double e_max) const
{
	// restrict integration to range where density can be nonzero
	double loge_min = std::max(minLog_, std::log(e_min));
	double loge_max = std::min(maxLog_, std::log(e_max));
	r_min = std::max(0., r_min);
	r_max = std::min(GetMaxRadius(), r_max);
	if (m > 1) {
		// Integrate dP/(dlogE dr^2) for numerical stability
		auto integrand = [this,d,ct,m](double r2, double loge)
		{
			double r = std::sqrt(r2);
			return std::exp(GetLog(d, ct, m, r, log_value(loge))
			    - std::log(2*r) + loge);
		};
		boost::array<double, 2> lo = {{r_min*r_min, loge_min}};
		boost::array<double, 2> hi = {{r_max*r_max, loge_max}};
		return I3MuonGun::Integrate(boost::function<double(double,double)>(integrand),
		    lo, hi, 1e-12, 1e-6, 10000);
	} else {
		// For single muons, dP/dr is a delta function at 0
		auto integrand = [this,d,ct,m](double loge)
		{
			return std::exp(GetLog(d, ct, m, 0, log_value(loge)) + loge);
		};
		return I3MuonGun::Integrate(boost::function<double(double)>(integrand), loge_min, loge_max);
	}
}

bool
SplineEnergyDistribution::operator==(const EnergyDistribution &o) const
{
	const SplineEnergyDistribution *other = dynamic_cast<const SplineEnergyDistribution*>(&o);
	if (!other)
		return false;
	else
		return (singles_ == other->singles_ && bundles_ == other->bundles_);
}

SplineEnergyDistribution::SplineEnergyDistribution(const std::string &singles, const std::string &bundles)
    : singles_(singles), bundles_(bundles)
{
	if (singles_.GetNDim() != 3u)
		log_fatal("'%s' does not appear to be a single-muon energy distribution", singles.c_str());
	if (bundles_.GetNDim() != 5u)
		log_fatal("'%s' does not appear to be a muon bundle energy distribution", bundles.c_str());
	SetMin(std::exp(std::max(singles_.GetExtents(2).first, bundles_.GetExtents(2).first)));
	SetMax(std::exp(std::min(singles_.GetExtents(2).second, bundles_.GetExtents(2).second)));
}

double
SplineEnergyDistribution::GetMaxRadius() const
{
	return bundles_.GetExtents(3).second;
}

double
SplineEnergyDistribution::GetLog(double depth, double cos_theta, 
    unsigned multiplicity, double radius, log_value log_energy) const
{
	double coords[5] = {cos_theta, depth, static_cast<double>(multiplicity),
	    radius, log_energy};
	double logprob;
	
	if (radius < 0 || radius > GetMaxRadius() ||
	    log_energy < minLog_ || log_energy > maxLog_) {
		return -std::numeric_limits<double>::infinity();
	} else if (multiplicity < 2) {
		coords[2] = coords[4];
		if (singles_.Eval(coords, &logprob) != 0)
			return -std::numeric_limits<double>::infinity();
	} else if (bundles_.Eval(coords, &logprob) != 0)
		return -std::numeric_limits<double>::infinity();
	
	// Bundle spline is fit to log(dP/dr^2 dlogE)
	if (multiplicity > 1)
		logprob += std::log(2*radius);
	
	return logprob;
}

std::vector<std::pair<double,double> >
SplineEnergyDistribution::Generate(I3RandomService &rng, double depth,
    double cos_theta, unsigned multiplicity, unsigned nsamples) const
{
	typedef double (Signature)(double, double);
	typedef EnsembleSampler<Signature> Sampler;
	
	// the number of walkers must be even, and at least twice the
	// dimensionality of the space
	const unsigned walkers = std::max(4u, nsamples + (nsamples % 2));
	std::vector<Sampler::array_type> initial_ensemble(walkers);
	{
		// Draw starting positions from the MUPAGE parameterization
		BMSSEnergyDistribution proposal;
		proposal.SetMin(GetMin());
		proposal.SetMax(GetMax());
		std::vector<std::pair<double, double> > samples =
		    proposal.Generate(rng, depth, cos_theta, multiplicity, walkers);
		for (unsigned i=0; i < walkers; i++) {
			initial_ensemble[i][0] = samples.at(i).first;
			initial_ensemble[i][1] = samples.at(i).second;
		}
	}
	
	auto log_posterior = [this,depth,cos_theta,multiplicity](double r, double e)
	{
		return this->GetLog(depth, cos_theta, multiplicity, r,
		    EnergyDistribution::log_value(std::log(e)));
	};
	Sampler sampler(log_posterior, initial_ensemble);
	
	// Run the sampler for a few cycles to make it independent of the initial
	// ensemble. Fewer than 50 or so burn-in steps is too small to reach the
	// stationary distribution, while more than 100 is a waste of time, as
	// measured with resources/test/test_sampling.py
	for (unsigned i=0; i < 64; i++)
		sampler.Sample(rng);
	
	// copy the current ensemble into the output
	std::vector<std::pair<double,double> > samples;
	samples.reserve(nsamples);
	const std::vector<Sampler::sample> &ensemble = sampler.Sample(rng);
	unsigned todo = std::min(walkers, nsamples - unsigned(samples.size()));
	for (unsigned j=0; j < todo; j++) {
		samples.push_back(std::make_pair(ensemble[j].point[0], ensemble[j].point[1]));
	}
	
	// check the acceptance rate for sanity.
	double acceptance_rate = sampler.GetAcceptanceRate();
	if (acceptance_rate < 0.2) {
		log_warn("Ensemble sampler accepted only %.0f%% of samples (too low). It may be stuck in a local maximum.", 100*acceptance_rate);
	} else if (acceptance_rate > 0.9) {
		log_warn("Ensemble sampler accepted %.0f%% of samples (too high). This is an expensive random walk.", 100*acceptance_rate);
	}
	
	return samples;
}

BMSSEnergyDistribution::BMSSEnergyDistribution() :
    beta_(0.42), g0_(-0.232), g1_(3.961), e0a_(0.0304), e0b_(0.359), e1a_(-0.0077), e1b_(0.659),
    a0_(0.0033), a1_(0.0079), b0a_(0.0407), b0b_(0.0283), b1a_(-0.312), b1b_(6.124),
    q0_(0.0543), q1_(-0.365), c0a_(-0.069), c0b_(0.488), c1_(-0.117),
    d0a_(-0.398), d0b_(3.955), d1a_(0.012), d1b_(-0.35)
{}

bool
BMSSEnergyDistribution::operator==(const EnergyDistribution &o) const
{
	return dynamic_cast<const BMSSEnergyDistribution*>(&o);
}

OffsetPowerLaw
BMSSEnergyDistribution::GetSpectrum(double depth, double cos_theta, unsigned m, double r) const
{
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
	double bX = beta_*h/cos_theta;
	double g, eps;
	if (m == 1) {
		g = g0_*log(h) + g1_;
		eps = (e0a_*exp(e0b_*h)/cos_theta + e1a_*h + e1b_)*I3Units::TeV;
	} else {
		m = std::min(m, 4u);
		double a = a0_*h + a1_;
		double b = (b0a_*m + b0b_)*h + (b1a_*m + b1b_);
		double q = q0_*h + q1_;
		g = a*r + b*(1 - 0.5*exp(q*r));
		double c = (c0a_*h + c0b_)*exp(c1_*r);
		double d = (d0a_*h + d0b_)*pow(r, d1a_*h + d1b_);
		eps = (c*acos(cos_theta) + d)*I3Units::TeV;
	}
	
	return OffsetPowerLaw(g, eps*(1-exp(-bX)), GetMin(), GetMax());
}

double
BMSSEnergyDistribution::GetLog(double depth, double cos_theta, 
    unsigned multiplicity, double radius, log_value log_energy) const
{
	
	return BMSSRadialDistribution().GetLog(depth, cos_theta, multiplicity, radius) +
	    GetSpectrum(depth, cos_theta, multiplicity, radius).GetLog(log_energy);
}

double
BMSSEnergyDistribution::GetMaxRadius() const
{
	return 250*I3Units::m;
}

std::vector<std::pair<double,double> >
BMSSEnergyDistribution::Generate(I3RandomService &rng, double depth,
    double cos_theta, unsigned multiplicity, unsigned samples) const
{
	std::vector<std::pair<double,double> > values;
	values.reserve(samples);
	std::pair<double, double> val;
	for (unsigned i=0; i < samples; i++) {
		val.first = BMSSRadialDistribution().Generate(rng, depth, cos_theta, multiplicity);
		val.second = GetSpectrum(depth, cos_theta, multiplicity, val.first).Generate(rng);
		values.push_back(val);
	}
	
	return values;
}


OffsetPowerLaw::OffsetPowerLaw() : gamma_(NAN), offset_(NAN), emin_(NAN), emax_(NAN)
{}

OffsetPowerLaw::OffsetPowerLaw(double gamma, double offset, double emin, double emax)
    : gamma_(gamma), offset_(offset), emin_(emin), emax_(emax)
{
	if (gamma <= 0)
		log_fatal("Power law index must be > 0");
	else if (gamma == 1) {
		nmin_ = std::log(emin + offset);
		nmax_ = std::log(emax + offset);
		norm_ = 1./(nmax_ - nmin_);
	} else {
		nmin_ = std::pow(emin + offset, 1-gamma);
		nmax_ = std::pow(emax + offset, 1-gamma);
		norm_ = (1-gamma)/(nmax_ - nmin_);
	}
	lognorm_ = std::log(norm_);
}

bool
OffsetPowerLaw::operator==(const OffsetPowerLaw &other) const
{
	return (gamma_ == other.gamma_ && offset_ == other.offset_
	    && emin_ == other.emin_ && emax_ == other.emax_ );
}

double
OffsetPowerLaw::operator()(double energy) const
{
	if (energy <= emax_ && energy >= emin_)
		return norm_*std::pow(energy + offset_, -gamma_);
	else
		return 0.;
}

double
OffsetPowerLaw::GetLog(double energy) const
{
	if (energy <= emax_ && energy >= emin_)
		return lognorm_ - gamma_*std::log(energy + offset_);
	else
		return -std::numeric_limits<double>::infinity();
}

double
OffsetPowerLaw::GetLog(EnergyDistribution::log_value log_energy) const
{
	return GetLog(std::exp(log_energy));
}

double
OffsetPowerLaw::Generate(I3RandomService &rng) const
{
	return InverseSurvivalFunction(rng.Uniform());
}

double
OffsetPowerLaw::InverseSurvivalFunction(double p) const
{
	if (gamma_ == 1)
		return std::exp((1-p)*(nmax_ - nmin_) + nmin_) - offset_;
	else
		return std::pow((1-p)*(nmax_ - nmin_) + nmin_, 1./(1.-gamma_)) - offset_;
}

template <typename Archive>
void
EnergyDistribution::serialize(Archive &ar __attribute__ ((unused)), unsigned version __attribute__ ((unused)))
{}
	
template <typename Archive>
void
SplineEnergyDistribution::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("EnergyDistribution", base_object<EnergyDistribution>(*this));
	ar & make_nvp("SingleEnergy", singles_);
	ar & make_nvp("BundleEnergy", bundles_);
}

template <typename Archive>
void
BMSSEnergyDistribution::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("EnergyDistribution", base_object<EnergyDistribution>(*this));
}

template <typename Archive>
void
OffsetPowerLaw::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Gamma", gamma_);
	ar & make_nvp("Offset", offset_);
	ar & make_nvp("MinEnergy", emin_);
	ar & make_nvp("MaxEnergy", emax_);
	
	*this = OffsetPowerLaw(gamma_, offset_, emin_, emax_);
}

}

I3_SERIALIZABLE(I3MuonGun::EnergyDistribution);
I3_SERIALIZABLE(I3MuonGun::SplineEnergyDistribution);
I3_SERIALIZABLE(I3MuonGun::BMSSEnergyDistribution);
I3_SERIALIZABLE(I3MuonGun::OffsetPowerLaw);
