/** $Id: StaticSurfaceInjector.cxx 136127 2015-08-12 13:54:32Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 136127 $
 * $Date: 2015-08-12 07:54:32 -0600 (Wed, 12 Aug 2015) $
 */

#include <MuonGun/I3MuonGun.h>
#include <MuonGun/StaticSurfaceInjector.h>
#include <MuonGun/Cylinder.h>
#include <dataclasses/I3Constants.h>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <icetray/I3Module.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <phys-services/I3RandomService.h>

namespace I3MuonGun {

namespace {

std::string
GetTablePath(const std::string &subpath)
{
	std::ostringstream tablePath;
	tablePath << getenv("I3_BUILD") << "/MuonGun/resources/tables/" << subpath;
	return tablePath.str();
}

}

template <typename Archive>
void
StaticSurfaceInjector::serialize(Archive &ar, unsigned version)
{
	if (version > 1)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Generator", base_object<Generator>(*this));
	if (version == 0) {
		CylinderPtr cyl;
		ar & make_nvp("Surface", cyl);
		surface_ = cyl;
	} else {
		ar & make_nvp("Surface", surface_);
	}
	ar & make_nvp("Flux", flux_);
	ar & make_nvp("EnergySpectrum", energyGenerator_);
	ar & make_nvp("RadialDistribution", radialDistribution_);
	ar & make_nvp("MaxFlux", maxFlux_);
	ar & make_nvp("TotalRate", totalRate_);
}

StaticSurfaceInjector::StaticSurfaceInjector()
{
	SetSurface(boost::make_shared<Cylinder>(1600, 800));
	
	FluxPtr flux = boost::make_shared<SplineFlux>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.single_flux.fits"),
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_flux.fits"));
	flux->SetMinMultiplicity(1);
	flux->SetMaxMultiplicity(1);
	SetFlux(flux);
	
	energyGenerator_ = boost::make_shared<OffsetPowerLaw>(2, 500., 50, 1e6);
	
	radialDistribution_ = boost::make_shared<SplineRadialDistribution>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.radius.fits"));
}

StaticSurfaceInjector::StaticSurfaceInjector(SamplingSurfacePtr surface, FluxPtr flux,
    boost::shared_ptr<OffsetPowerLaw> edist, RadialDistributionPtr rdist)
{
	SetSurface(surface);
	SetFlux(flux);
	energyGenerator_ = edist;
	radialDistribution_ = rdist;
}

GenerationProbabilityPtr
StaticSurfaceInjector::Clone() const
{
	return boost::make_shared<StaticSurfaceInjector>(*this);
}

bool
StaticSurfaceInjector::IsCompatible(GenerationProbabilityConstPtr o) const
{
	boost::shared_ptr<const StaticSurfaceInjector> other
	    = boost::dynamic_pointer_cast<const StaticSurfaceInjector>(o);
	if (!other)
		return false;
	return (*surface_ == *(other->surface_) && *flux_ == *(other->flux_)
	    && *radialDistribution_ == *(other->radialDistribution_)
	    && *energyGenerator_ == *(other->energyGenerator_));
}

void
StaticSurfaceInjector::SetSurface(SamplingSurfacePtr p)
{
	surface_ = p;
	totalRate_ = NAN;
	zenithNorm_ = NAN;
	CalculateMaxFlux();
}

void
StaticSurfaceInjector::SetFlux(FluxPtr p)
{
	flux_ = p;
	totalRate_ = NAN;
	zenithNorm_ = NAN;
	CalculateMaxFlux();
}

void
StaticSurfaceInjector::CalculateMaxFlux()
{
	if (surface_ && flux_)
		maxFlux_ = (*flux_)(surface_->GetMinDepth(), 1., 1u)*surface_->GetMaximumArea();
}

double
StaticSurfaceInjector::GetTotalRate() const
{
	if (std::isnan(totalRate_) && surface_ && flux_) {
		totalRate_ = 0;
		for (unsigned m = flux_->GetMinMultiplicity(); m <= flux_->GetMaxMultiplicity(); m++) {
			totalRate_ += surface_->IntegrateFlux(boost::bind(boost::cref(*flux_), _1, _2, m));
		}
	}
	return totalRate_;
}

double
StaticSurfaceInjector::GetZenithNorm() const
{
	if (std::isnan(zenithNorm_) && surface_ && flux_) {
		zenithNorm_ = 0;
		for (unsigned m = flux_->GetMinMultiplicity(); m <= flux_->GetMaxMultiplicity(); m++) {
			zenithNorm_ += Integrate(boost::bind(boost::cref(*flux_), surface_->GetMinDepth(), _1, m), 0, 1);
		}
		zenithNorm_ = std::log(zenithNorm_);
	}
	return zenithNorm_;
}

void
StaticSurfaceInjector::GenerateAxis(I3RandomService &rng, std::pair<I3Particle, unsigned> &axis) const
{
	// Choose a direction and impact position from a uniform flux through
	// the sampling surface, then pick a multiplicity. Accept zenith angles
	// and multiplicities at a rate proportional to the flux at the shallowest depth.
	I3Direction dir;
	I3Position pos;
	unsigned m;
	double flux;
	// For any sane distribution the maximum flux is for single vertical muons
	// at minimum depth
	double maxflux = (*flux_)(surface_->GetMinDepth(), 1., flux_->GetMinMultiplicity());
	do {
		surface_->SampleImpactRay(pos, dir, rng);
		m = rng.Integer(flux_->GetMaxMultiplicity() - flux_->GetMinMultiplicity())
		    + flux_->GetMinMultiplicity();
		flux = (*flux_)(surface_->GetMinDepth(), cos(dir.GetZenith()), m);
	} while (rng.Uniform(0., maxflux) > flux);
	
	axis.first.SetPos(pos);
	axis.first.SetDir(dir);
	axis.first.SetShape(I3Particle::Primary);
	axis.first.SetLocationType(I3Particle::Anywhere);
	axis.first.SetType(I3Particle::unknown);
	axis.first.SetTime(0.);
	axis.second = m;
}

void
StaticSurfaceInjector::FillMCTree(I3RandomService &rng,
    const std::pair<I3Particle, unsigned> &axis,
    I3MCTree &mctree, BundleConfiguration &bundlespec) const
{
	const I3Particle &primary = axis.first;
	I3MCTreeUtils::AddPrimary(mctree, primary);
	double h = GetDepth(primary.GetPos().GetZ());
	double coszen = cos(primary.GetDir().GetZenith());
	
	unsigned m = axis.second;
	for (unsigned i=0; i < m; i++) {
		double radius = 0., azimuth = 0.;
		if (m > 1u) {
			radius = radialDistribution_->Generate(rng, h, coszen, m);
			azimuth = rng.Uniform(0., 2*M_PI);
		}
		
		I3Particle track = CreateParallelTrack(radius, azimuth, *surface_, primary);
		
		track.SetEnergy(energyGenerator_->Generate(rng));
		I3MCTreeUtils::AppendChild(mctree, primary, track);
		bundlespec.push_back(BundleEntry(radius, track.GetEnergy()));
	}
}

void
StaticSurfaceInjector::Generate(I3RandomService &rng, I3MCTree &mctree, BundleConfiguration &bundlespec) const
{
	std::pair<I3Particle, unsigned> axis;
	GenerateAxis(rng, axis);
	FillMCTree(rng, axis, mctree, bundlespec);
}

double
StaticSurfaceInjector::GetLogGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundlespec) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
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
	
	return logprob - std::log(surface_->GetAcceptance());
}

}

I3_SERIALIZABLE(I3MuonGun::StaticSurfaceInjector);
