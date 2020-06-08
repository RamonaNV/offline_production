/** $Id: NaturalRateInjector.cxx 165001 2018-08-27 09:22:30Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 165001 $
 * $Date: 2018-08-27 03:22:30 -0600 (Mon, 27 Aug 2018) $
 */

#include <MuonGun/I3MuonGun.h>
#include <MuonGun/NaturalRateInjector.h>
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
NaturalRateInjector::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Generator", base_object<Generator>(*this));
	ar & make_nvp("Surface", surface_);
	
	ar & make_nvp("Flux", flux_);
	ar & make_nvp("EnergyRadiusDistribution", energyDistribution_);
	ar & make_nvp("TotalRate", totalRate_);
}

NaturalRateInjector::NaturalRateInjector()
{
	SetSurface(boost::make_shared<Cylinder>(1600, 800));
	
	SetFlux(boost::make_shared<SplineFlux>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.single_flux.fits"),
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_flux.fits")));
	
	SetEnergyDistribution(boost::make_shared<SplineEnergyDistribution>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.single_energy.fits"),
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_energy.fits")));
}

NaturalRateInjector::NaturalRateInjector(SamplingSurfacePtr surface, FluxPtr flux,
    EnergyDistributionPtr edist)
{
	SetSurface(surface);
	SetFlux(flux);
	SetEnergyDistribution(edist);
}

GenerationProbabilityPtr
NaturalRateInjector::Clone() const
{
	return boost::make_shared<NaturalRateInjector>(*this);
}

bool
NaturalRateInjector::IsCompatible(GenerationProbabilityConstPtr o) const
{
	boost::shared_ptr<const NaturalRateInjector> other
	    = boost::dynamic_pointer_cast<const NaturalRateInjector>(o);
	if (!other)
		return false;
	return (*surface_ == *(other->surface_) && *flux_ == *(other->flux_)
	    && *energyDistribution_ == *(other->energyDistribution_));
}

void
NaturalRateInjector::SetSurface(SamplingSurfacePtr p)
{
	assert(p);
	surface_ = p;
	totalRate_ = NAN;
}

void
NaturalRateInjector::SetFlux(FluxPtr p)
{
	assert(p);
	flux_ = p;
	totalRate_ = NAN;
}

void
NaturalRateInjector::SetEnergyDistribution(EnergyDistributionPtr e)
{
	assert(e);
	
	double norm = e->Integrate(surface_->GetMinDepth(), 1., flux_->GetMinMultiplicity(), 0, 300, e->GetMin(), e->GetMax());
	if (std::abs(norm-1) > 1e-1)
		log_fatal_stream("The provided energy distribution is not normalized "
		    "(integrates to "<<norm<<"). This could be because the minimum and "
		    "maximum energies were changed from their defaults.");
	energyDistribution_ = e;
}


double
NaturalRateInjector::GetTotalRate() const
{
	if (std::isnan(totalRate_) && surface_ && flux_) {
		totalRate_ = 0;
		for (unsigned m = flux_->GetMinMultiplicity(); m <= flux_->GetMaxMultiplicity(); m++) {
			totalRate_ += surface_->IntegrateFlux(boost::bind(boost::cref(*flux_), _1, _2, m));
		}
	}
	return totalRate_;
}

void
NaturalRateInjector::GenerateAxis(I3RandomService &rng, std::pair<I3Particle, unsigned> &axis) const
{
	// Choose a direction and impact position from a uniform flux through
	// the sampling surface, then pick a multiplicity. Accept zenith angles
	// and multiplicities at a rate proportional to the flux at the chosen
	// depth.
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
		flux = (*flux_)(GetDepth(pos.GetZ()), cos(dir.GetZenith()), m);
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
NaturalRateInjector::FillMCTree(I3RandomService &rng,
    const std::pair<I3Particle, unsigned> &axis,
    I3MCTree &mctree, BundleConfiguration &bundlespec) const
{
	const I3Particle &primary = axis.first;
	I3MCTreeUtils::AddPrimary(mctree, primary);
	double h = GetDepth(primary.GetPos().GetZ());
	double coszen = cos(primary.GetDir().GetZenith());
	
	unsigned m = axis.second;
	auto tracks = energyDistribution_->Generate(rng, h, coszen, m, m);
	for (auto &radius_energy : tracks) {
		I3Particle track = CreateParallelTrack(radius_energy.first, rng.Uniform(0., 2*M_PI), *surface_, primary);
		track.SetEnergy(radius_energy.second);
		I3MCTreeUtils::AppendChild(mctree, primary, track);
		bundlespec.push_back(BundleEntry(radius_energy.first, radius_energy.second));
	}
}

void
NaturalRateInjector::Generate(I3RandomService &rng, I3MCTree &mctree, BundleConfiguration &bundlespec) const
{
	std::pair<I3Particle, unsigned> axis;
	GenerateAxis(rng, axis);
	FillMCTree(rng, axis, mctree, bundlespec);
}

double
NaturalRateInjector::GetLogGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundlespec) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return -std::numeric_limits<double>::infinity();
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = static_cast<unsigned>(bundlespec.size());

	// We used the flux to do rejection sampling in depth, zenith, and
	// multiplicity. Evaluate the properly-normalized PDF here.
	double logprob = flux_->GetLog(h, coszen, m) - std::log(GetTotalRate());
	BOOST_FOREACH(const BundleEntry &track, bundlespec) {
		logprob += energyDistribution_->GetLog(h, coszen, m, track.radius,
		    EnergyDistribution::log_value(std::log(track.energy)));
	}
	
	return logprob;
}

}

I3_SERIALIZABLE(I3MuonGun::NaturalRateInjector);
