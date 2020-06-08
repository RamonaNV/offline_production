/** $Id: RadialDistribution.cxx 159352 2017-11-07 22:43:41Z juancarlos $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 159352 $
 * $Date: 2017-11-07 15:43:41 -0700 (Tue, 07 Nov 2017) $
 */

#include <MuonGun/RadialDistribution.h>
#include <phys-services/I3RandomService.h>
#include <icetray/I3Units.h>

namespace I3MuonGun {

RadialDistribution::~RadialDistribution() {}

double
RadialDistribution::operator()(double depth, double cos_theta,
    unsigned N, double radius) const
{
	if (N < 2)
		return 1.;
	else
		return std::exp(GetLog(depth, cos_theta, N, radius));
}

BMSSRadialDistribution::BMSSRadialDistribution() : rho0a_(-1.786), rho0b_(28.26),
    rho1_(-1.06), theta0_(1.3), f_(10.4), alpha0a_(-0.448), alpha0b_(4.969),
    alpha1a_(0.0194), alpha1b_(0.276), rmax_(250*I3Units::m) {};

double
BMSSRadialDistribution::GetMeanRadius(double h, double theta, unsigned N) const
{
	unsigned n = std::min(N, 4u);
	return ((rho0a_*n + rho0b_)*pow(h,rho1_))/(exp((theta-theta0_)*f_)+1.);
}

double
BMSSRadialDistribution::GetShapeParameter(double h, double theta __attribute__ ((unused)), unsigned N) const
{
	unsigned n = std::min(N, 4u);
	return (alpha0a_*n + alpha0b_)*exp(h*(alpha1a_*n + alpha1b_));
}

double
BMSSRadialDistribution::GetGenerationProbability(double R, double a, double radius) const
{
	double R0 = R*(a-3)/2.;
	
	return (a-1)*(a-2)*pow(R0, a-2)*(radius/pow(radius+R0, a));
}

double
BMSSRadialDistribution::GetLog(double depth, double cos_theta,
    unsigned N, double radius) const
{
	if (!(N > 1))
		return 0.;
	
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
	double theta = acos(cos_theta);
	
	return std::log(GetGenerationProbability(GetMeanRadius(h, theta, N), GetShapeParameter(h, theta, N), radius));
}

double
BMSSRadialDistribution::Generate(I3RandomService &rng, double depth, double cos_theta,
    unsigned N) const
{
	if (!(N > 1))
		return 0.;
	
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
	double theta = acos(cos_theta);
	double R = GetMeanRadius(h, theta, N);
	double a = GetShapeParameter(h, theta, N);
	
	double peak_radius = std::max(0., R*(a-3)/(2.*(a-1)));
	if (!std::isfinite(peak_radius))
		log_fatal("Peak radius is not finite!");
	double max_prob = GetGenerationProbability(R, a, peak_radius);
	if (!std::isfinite(max_prob))
		log_fatal("Peak probability is not finite!");
	double r;
	do {
		r = rng.Uniform(rmax_);
	} while (rng.Uniform(max_prob) > GetGenerationProbability(R, a, r));

	return r;
}

bool
BMSSRadialDistribution::operator==(const RadialDistribution &o) const
{
	return dynamic_cast<const BMSSRadialDistribution*>(&o);
}

SplineRadialDistribution::SplineRadialDistribution(const std::string &path)
    : spline_(path) {}

double
SplineRadialDistribution::GetLog(double depth, double cos_theta,
    unsigned N, double radius) const
{
	double coords[4] = {cos_theta, depth, static_cast<double>(N), radius};
	double logprob;
	
	if (spline_.Eval(coords, &logprob) != 0)
		return -std::numeric_limits<double>::infinity();
	else
		// Spline is fit to log(dP/dr^2)
		return std::log(2*radius) + logprob;
}

double
SplineRadialDistribution::Generate(I3RandomService &rng, double depth,
    double cos_theta, unsigned N) const
{
	double radius, logprob, maxprob; (void) radius;
	std::pair<double, double> extent = spline_.GetExtents(3);
	double coords[4] = {cos_theta, depth, static_cast<double>(N), extent.first};
	if (spline_.Eval(coords, &maxprob) != 0)
		maxprob = -std::numeric_limits<double>::infinity();
	
	// The spline is fit to log(dP/dr^2) as a function of r,
	// so we generate proposals uniformly in r^2, then take
	// a square root to evaluate.
	do {
		coords[3] = std::sqrt(rng.Uniform(extent.first*extent.first,
		    extent.second*extent.second));
		if (spline_.Eval(coords, &logprob) != 0)
			logprob = -std::numeric_limits<double>::infinity();
	} while (std::log(rng.Uniform()) > logprob - maxprob);
	
	return coords[3];
}

bool
SplineRadialDistribution::operator==(const RadialDistribution &o) const
{
	const SplineRadialDistribution *other = dynamic_cast<const SplineRadialDistribution*>(&o);
	if (!other)
		return false;
	else
		return (spline_ == other->spline_);
}

template <typename Archive>
void
RadialDistribution::serialize(Archive &ar __attribute__((unused)), unsigned __attribute__((unused)))
{}

template <typename Archive>
void
SplineRadialDistribution::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("RadialDistribution", base_object<RadialDistribution>(*this));
	ar & make_nvp("SplineTable", spline_);
}

template <typename Archive>
void
BMSSRadialDistribution::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("RadialDistribution", base_object<RadialDistribution>(*this));
}

}

I3_SERIALIZABLE(I3MuonGun::RadialDistribution);
I3_SERIALIZABLE(I3MuonGun::SplineRadialDistribution);
I3_SERIALIZABLE(I3MuonGun::BMSSRadialDistribution);

