
#include <MuonGun/Floodlight.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <MuonGun/Cylinder.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

#include <boost/foreach.hpp>

namespace I3MuonGun {

template <typename Archive>
void
Floodlight::serialize(Archive &ar, unsigned version)
{
	if (version > 1)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Generator", base_object<Generator>(*this));
	ar & make_nvp("Surface", surface_);
	ar & make_nvp("EnergySpectrum", energyGenerator_);
	if (version == 0) {
		zenith_range_ = std::make_pair(-1., 1.);
	} else {
		ar & make_nvp("CosZenithRange", zenith_range_);
	}
	log_acceptance_ = std::log(surface_->GetAcceptance(zenith_range_.first, zenith_range_.second));
	
}

Floodlight::Floodlight(SamplingSurfacePtr surface, boost::shared_ptr<OffsetPowerLaw> energyGenerator,
    double cosMin, double cosMax) : surface_(surface), energyGenerator_(energyGenerator),
    zenith_range_(cosMin, cosMax)
{
	if (!surface_)
		surface_ = boost::make_shared<Cylinder>(1000, 600, I3Position(31.25, 19.64, 0));
	if (!energyGenerator_)
		energyGenerator_ = boost::make_shared<OffsetPowerLaw>(1, 0., 5e2, 1e7);
	i3_assert(zenith_range_.second > zenith_range_.first);
	i3_assert(std::abs(zenith_range_.first) <= 1);
	i3_assert(std::abs(zenith_range_.second) <= 1);
	log_acceptance_ = std::log(surface_->GetAcceptance(zenith_range_.first, zenith_range_.second));
	i3_assert(std::isfinite(log_acceptance_));
}

GenerationProbabilityPtr
Floodlight::Clone() const
{
	return GenerationProbabilityPtr(new Floodlight(*this));
}

bool
Floodlight::IsCompatible(GenerationProbabilityConstPtr o) const
{
	boost::shared_ptr<const Floodlight> other = boost::dynamic_pointer_cast<const Floodlight>(o);
	if (!other)
		return false;
	else
		return (*surface_ == *(other->surface_)
		    && zenith_range_ == other->zenith_range_
		    && *energyGenerator_ == *(other->energyGenerator_));
}

void
Floodlight::Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle __attribute__((unused))) const
{
	I3Direction dir;
	I3Position pos;
	surface_->SampleImpactRay(pos, dir, rng, zenith_range_.first, zenith_range_.second);
	
	I3Particle primary;
	primary.SetDir(dir);
	primary.SetPos(pos);
	primary.SetTime(0);
	primary.SetEnergy(energyGenerator_->Generate(rng));
	
	I3Particle muon = primary.Clone();
	muon.SetType(I3Particle::MuMinus);
	muon.SetLocationType(I3Particle::InIce);
	muon.SetShape(I3Particle::Null);
	
	I3MCTreeUtils::AddPrimary(tree, primary);
	I3MCTreeUtils::AppendChild(tree, primary, muon);
}

double
Floodlight::GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// Bail if the axis doesn't intersect the surface, or there's more than 1 muon.
	double ct = std::cos(axis.GetDir().GetZenith());
	if (!std::isfinite(steps.first) || bundle.size() != 1
	    || ct < zenith_range_.first
	    || ct > zenith_range_.second)
		return -std::numeric_limits<double>::infinity();
	
	return energyGenerator_->GetLog(bundle.front().energy) - log_acceptance_;
}

}

I3_SERIALIZABLE(I3MuonGun::Floodlight);
