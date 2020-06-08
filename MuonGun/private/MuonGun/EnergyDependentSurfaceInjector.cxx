/** $Id: EnergyDependentSurfaceInjector.cxx 136127 2015-08-12 13:54:32Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 136127 $
 * $Date: 2015-08-12 07:54:32 -0600 (Wed, 12 Aug 2015) $
 */

#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/SamplingSurface.h>
#include <MuonGun/Cylinder.h>
#include <MuonGun/Flux.h>
#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/RadialDistribution.h>
#include <dataclasses/I3Constants.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <phys-services/I3RandomService.h>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

namespace I3MuonGun {

EnergyDependentSurfaceInjector::EnergyDependentSurfaceInjector(CylinderPtr surface, FluxPtr flux,
    boost::shared_ptr<OffsetPowerLaw> energies, RadialDistributionPtr radius, SurfaceScalingFunctionPtr scaling)
    : StaticSurfaceInjector(surface, flux, energies, radius), scalingFunction_(scaling)
{
}

GenerationProbabilityPtr
EnergyDependentSurfaceInjector::Clone() const
{
	return boost::make_shared<EnergyDependentSurfaceInjector>(*this);
}

bool
EnergyDependentSurfaceInjector::IsCompatible(GenerationProbabilityConstPtr o) const
{
	if (!StaticSurfaceInjector::IsCompatible(o))
		return false;
	boost::shared_ptr<const EnergyDependentSurfaceInjector> other
	    = boost::dynamic_pointer_cast<const EnergyDependentSurfaceInjector>(o);
	if (!other)
		return false;
	else
		return *scalingFunction_ == *(other->scalingFunction_);
}

void
EnergyDependentSurfaceInjector::Generate(I3RandomService &rng, I3MCTree &tree,
    BundleConfiguration &bundle) const
{
	I3Direction dir;
	I3Position pos;
	unsigned m;
	double flux;
	double maxflux = (*flux_)(surface_->GetMinDepth(), 1., flux_->GetMinMultiplicity());
	double h, coszen;
	SamplingSurfaceConstPtr surface;
	do {
		// Choose a multiplicity
		bundle.clear();
		m = rng.Integer(flux_->GetMaxMultiplicity() - flux_->GetMinMultiplicity())
		    + flux_->GetMinMultiplicity();
		// Choose an ensemble of energies
		for (unsigned i=0; i < m; i++)
			bundle.push_back(BundleEntry(0., energyGenerator_->Generate(rng)));
		bundle.sort();
		
		// Choose target surface based on highest-energy muon
		surface = GetTargetSurface(bundle.front().energy);
		// Sample an impact point on the target surface
		surface->SampleImpactRay(pos, dir, rng);
		h = GetDepth(pos.GetZ());
		coszen = cos(dir.GetZenith());
		
		// Snap the impact point back to the injection surface
		std::pair<double, double> steps = surface_->GetIntersection(pos, dir);
		if (!(steps.first <= 0))
			log_fatal("The target point is outside the injection surface!");
		pos.SetX(pos.GetX() + steps.first*dir.GetX());
		pos.SetY(pos.GetY() + steps.first*dir.GetY());
		pos.SetZ(pos.GetZ() + steps.first*dir.GetZ());
		
		// Accept or reject the chosen zenith angle and multiplicity
		flux = (*flux_)(surface_->GetMinDepth(), cos(dir.GetZenith()), m);
	} while (rng.Uniform(0., maxflux) > flux);

	I3Particle primary;
	primary.SetPos(pos);
	primary.SetDir(dir);
	primary.SetShape(I3Particle::Primary);
	primary.SetLocationType(I3Particle::Anywhere);
	primary.SetType(I3Particle::unknown);
	primary.SetTime(0.);
	I3MCTreeUtils::AddPrimary(tree, primary);
	
	// For each muon, draw a radial offset and add an entry
	// to the MCTree
	BOOST_FOREACH(BundleConfiguration::value_type &bspec, bundle) {
		double radius = 0., azimuth = 0.;
		if (m > 1u) {
			radius = radialDistribution_->Generate(rng, h, coszen, m);
			azimuth = rng.Uniform(0., 2*M_PI);
		}
		
		I3Particle track = CreateParallelTrack(radius, azimuth, *surface, primary);
		track.SetEnergy(bspec.energy);
		bspec.radius = radius;
		I3MCTreeUtils::AppendChild(tree, primary, track);
	}
}

SamplingSurfacePtr
EnergyDependentSurfaceInjector::GetTargetSurface(double energy) const
{
	if (scalingFunction_)
		return scalingFunction_->GetSurface(energy);
	else
		return boost::make_shared<Cylinder>(1600, 800);
}

double
EnergyDependentSurfaceInjector::GetTotalRate(SamplingSurfaceConstPtr surface) const
{
	double total_rate = 0.;
	for (unsigned m = flux_->GetMinMultiplicity(); m <= flux_->GetMaxMultiplicity(); m++)
		total_rate += surface->IntegrateFlux(boost::bind(boost::cref(*flux_), _1, _2, m));
	return total_rate;
}

double
EnergyDependentSurfaceInjector::GetLogGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundlespec) const
{
	// Entries are sorted in descending order of energy, so the
	// "minimum" entry has the maximum energy
	SamplingSurfaceConstPtr surface =
	    GetTargetSurface(std::min_element(bundlespec.begin(), bundlespec.end())->energy);
	std::pair<double, double> steps =
	    surface->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the target surface. Bail.
	if (!std::isfinite(steps.first))
		return -std::numeric_limits<double>::infinity();
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = static_cast<unsigned>(bundlespec.size());

	// We used the flux to do rejection sampling in zenith and multiplicity. Evaluate
	// the properly-normalized PDF here.
	double logprob = flux_->GetLog(surface_->GetMinDepth(), coszen, m) - GetZenithNorm();
	BOOST_FOREACH(const BundleEntry &track, bundlespec) {
		if (m > 1)
			logprob += radialDistribution_->GetLog(h, coszen, m, track.radius);
		logprob += energyGenerator_->GetLog(track.energy);
	}
	
	// We only distributed events over the target surface, not the entire injection surface
	return logprob - std::log(surface->GetAcceptance());
}

template <typename Archive>
void
EnergyDependentSurfaceInjector::serialize(Archive &ar, unsigned)
{
	ar & make_nvp("StaticSurfaceInjector", base_object<StaticSurfaceInjector>(*this));
	ar & make_nvp("ScalingFunction", scalingFunction_);
}

SurfaceScalingFunction::~SurfaceScalingFunction() {}

template <typename Archive>
void
SurfaceScalingFunction::serialize(Archive &ar __attribute__((unused)), unsigned version __attribute__((unused)))
{}

ConstantSurfaceScalingFunction::ConstantSurfaceScalingFunction() {}
ConstantSurfaceScalingFunction::ConstantSurfaceScalingFunction(SamplingSurfacePtr surface) : surface_(surface) {}
ConstantSurfaceScalingFunction::~ConstantSurfaceScalingFunction() {}

template <typename Archive>
void
ConstantSurfaceScalingFunction::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("SurfaceScalingFunction", base_object<SurfaceScalingFunction>(*this));
	ar & make_nvp("Surface", surface_);
}

SamplingSurfacePtr
ConstantSurfaceScalingFunction::GetSurface(double energy __attribute__((unused))) const { return surface_; }

bool
ConstantSurfaceScalingFunction::operator==(const SurfaceScalingFunction &o) const
{
	const ConstantSurfaceScalingFunction *other = dynamic_cast<const ConstantSurfaceScalingFunction*>(&o);
	if (!other)
		return false;
	else
		return (*surface_ == *(other->surface_));
}

template <typename Archive>
void
BasicSurfaceScalingFunction::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("SurfaceScalingFunction", base_object<SurfaceScalingFunction>(*this));
	ar & make_nvp("Scale", scale_);
	ar & make_nvp("EnergyScale", energyScale_);
	ar & make_nvp("Offset", offset_);
	ar & make_nvp("Power", power_);
	ar & make_nvp("RadiusBounds", rBounds_);
	ar & make_nvp("ZBounds", zBounds_);
	ar & make_nvp("CenterBounds", centerBounds_);
}

BasicSurfaceScalingFunction::BasicSurfaceScalingFunction() :
    scale_(800., 240957.5), energyScale_(4., 4.), offset_(3.778, 3.622), power_(1.10, 2.23),
    rBounds_(0., 525.), zBounds_(-500., 400.),
    centerBounds_(pair(46.29,-34.88), pair(31.25, 19.64))
{}

BasicSurfaceScalingFunction::~BasicSurfaceScalingFunction() {}

double
BasicSurfaceScalingFunction::GetMargin(double logenergy, double scale, double offset, double power) const
{
	if (logenergy < offset)
		return pow(scale*(offset - logenergy), 1./power);
	else
		return 0.;
}

SamplingSurfacePtr
BasicSurfaceScalingFunction::GetSurface(double energy) const
{
	// Shrink the cylinder down by an energy-dependent amount
	double z = std::max(zBounds_.second -
	    GetMargin(std::log10(energy/energyScale_.first), scale_.first, offset_.first, power_.first), zBounds_.first);
	
	// Shrink the sides of the cylinder
	double r = std::max(rBounds_.second -
	    GetMargin(std::log10(energy/energyScale_.second), scale_.second, offset_.second, power_.second), rBounds_.first);
	
	// Move the center smoothly between the configured bounds
	double hscale = r/(rBounds_.second-rBounds_.first);
	I3Position center(centerBounds_.first.first + hscale*(centerBounds_.second.first-centerBounds_.first.first),
	    centerBounds_.first.second + hscale*(centerBounds_.second.second-centerBounds_.first.second),
	    (zBounds_.first + z)/2.);
	
	return boost::make_shared<Cylinder>(z-zBounds_.first, r, center);
}

void
BasicSurfaceScalingFunction::SetCapScaling(double energyScale, double scale, double offset, double power)
{
	energyScale_.first = energyScale;
	scale_.first = scale;
	offset_.first = offset;
	power_.first = power;
}

void
BasicSurfaceScalingFunction::SetSideScaling(double energyScale, double scale, double offset, double power)
{
	energyScale_.second = energyScale;
	scale_.second = scale;
	offset_.second = offset;
	power_.second = power;
}

void
BasicSurfaceScalingFunction::SetRadiusBounds(double rmin, double rmax)
{
	rBounds_ = std::make_pair(rmin, rmax);
}

void
BasicSurfaceScalingFunction::SetZBounds(double zmin, double zmax)
{
	zBounds_ = std::make_pair(zmin, zmax);
}

void
BasicSurfaceScalingFunction::SetCenterMin(double x, double y)
{
	centerBounds_.first = std::make_pair(x, y);
}

void
BasicSurfaceScalingFunction::SetCenterMax(double x, double y)
{
	centerBounds_.second = std::make_pair(x, y);
}

bool
BasicSurfaceScalingFunction::operator==(const SurfaceScalingFunction &o) const
{
	const BasicSurfaceScalingFunction *other = dynamic_cast<const BasicSurfaceScalingFunction*>(&o);
	if (!other)
		return false;
	else
		return (scale_ == other->scale_ &&
		    energyScale_ == other->energyScale_ &&
		    offset_ == other->offset_ &&
		    power_ == other->power_ &&
		    rBounds_ == other->rBounds_ &&
		    zBounds_ == other->zBounds_ &&
		    centerBounds_ == other->centerBounds_);
}

}

I3_SERIALIZABLE(I3MuonGun::SurfaceScalingFunction);
I3_SERIALIZABLE(I3MuonGun::ConstantSurfaceScalingFunction);
I3_SERIALIZABLE(I3MuonGun::BasicSurfaceScalingFunction);
I3_SERIALIZABLE(I3MuonGun::EnergyDependentSurfaceInjector);
