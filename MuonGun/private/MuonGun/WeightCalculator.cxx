/** $Id: WeightCalculator.cxx 178081 2019-12-18 19:13:26Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 178081 $
 * $Date: 2019-12-18 12:13:26 -0700 (Wed, 18 Dec 2019) $
 */

#include <MuonGun/WeightCalculator.h>
#include <MuonGun/Generator.h>
#include <MuonGun/SamplingSurface.h>
#include <MuonGun/Cylinder.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/Track.h>
#include <boost/foreach.hpp>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
#include <simclasses/I3MMCTrack.h>
#include <phys-services/I3Calculator.h>
#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace I3MuonGun {

double
WeightCalculator::GetWeight(const I3Particle &axis, const BundleConfiguration &bundlespec) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return 0.;
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = unsigned(bundlespec.size());
	
	double rate = flux_->GetLog(h, coszen, m) - generator_->GetLogGeneratedEvents(axis, bundlespec);
	
	BOOST_FOREACH(const BundleEntry &track, bundlespec){
		double logprob = energy_->GetLog(h, coszen, m, track.radius,
		    EnergyDistribution::log_value(std::log(track.energy)));
		if (!std::isfinite(logprob)){
                    log_warn("Log Energy weight of a least one muon is -inf, weight will be 0!");
                }
		rate += logprob;
        }
	// assert(std::isfinite(std::exp(rate)));
	return std::exp(rate);
}

// Possibly throw-away utility function: "track" muons to a fixed surface using the
// same method as WeightCalculatorModule
std::vector<I3Particle>
GetMuonsAtSurface(I3FramePtr frame, I3Surfaces::SurfaceConstPtr surface)
{
	std::vector<I3Particle> final_states;
	
	I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
	I3MMCTrackListConstPtr mmctracks = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");
	if (!mctree)
		log_fatal("I3MCTree missing!");
	if (!mmctracks)
		log_fatal("I3MMCTrackList missing!");
	BOOST_FOREACH(const Track &track, Track::Harvest(*mctree, *mmctracks)) {
		std::pair<double, double> steps =
		    surface->GetIntersection(track.GetPos(), track.GetDir());
		double energy = track.GetEnergy(steps.first);
		if (steps.first >= 0 && energy > 0) {
			final_states.push_back(track);
			I3Particle &p = final_states.back();
			p.SetEnergy(energy);
			p.SetPos(track.GetPos(steps.first));
			p.SetTime(track.GetTime(steps.first));
			p.SetLength(track.GetLength()-steps.first);
		}
	}
	
	return final_states;
}

namespace {

namespace ublas = boost::numeric::ublas;
typedef ublas::bounded_vector<double, 3> vector;

vector
make_vector(double x, double y, double z)
{
	vector v(3);
	v[0] = x; v[1] = y; v[2] = z;
	
	return v;
}

vector
make_vector(const I3Direction &dir)
{
	return make_vector(dir.GetX(), dir.GetY(), dir.GetZ());
}

inline vector
subtract(const I3Position &p1, const I3Position &p2)
{
	return make_vector(p1.GetX()-p2.GetX(), p1.GetY()-p2.GetY(), p1.GetZ()-p2.GetZ());
}

inline double
GetRadius(const I3Particle &axis, const I3Position &pos)
{
        vector r = subtract(pos,axis.GetPos());
	double l = ublas::inner_prod(make_vector(axis.GetDir()), r);
	
	return sqrt(std::max(0., ublas::inner_prod(r, r) - l*l));
}

}

MuonBundleConverter::MuonBundleConverter(size_t maxMultiplicity, SamplingSurfaceConstPtr surface)
    : maxMultiplicity_(maxMultiplicity),
    surface_(surface ? surface : boost::make_shared<Cylinder>(1600, 800))
{}

I3TableRowDescriptionPtr
MuonBundleConverter::CreateDescription(const I3MCTree&)
{
	I3TableRowDescriptionPtr desc(new I3TableRowDescription());
	
	desc->AddField<uint32_t>("multiplicity", "", "Number of muons in the bundle");
	desc->AddField<float>("depth", "km", "Vertical depth of intersection with the sampling surface");
	desc->AddField<float>("cos_theta", "", "Cosine of the shower zenith angle");
	desc->AddField<float>("energy", "GeV", "Muon energy at sampling surface",
	    maxMultiplicity_);
	desc->AddField<float>("radius", "m", "Perpendicular distance from of track "
	    "from the bundle axis at the sampling surface", maxMultiplicity_);
	
	return desc;
}

size_t
MuonBundleConverter::FillRows(const I3MCTree &mctree, I3TableRowPtr rows)
{
	I3MMCTrackListConstPtr mmctracks = currentFrame_->Get<I3MMCTrackListConstPtr>("MMCTrackList");
	if (!mmctracks)
		log_fatal("I3MMCTrackList missing!");
	
	const I3MCTree::const_iterator primary = mctree.begin();
	std::pair<double, double> primary_steps =
	    surface_->GetIntersection(primary->GetPos(), primary->GetDir());
	if (primary_steps.first > 0) {
		rows->Set<float>("depth", float(GetDepth(primary->GetPos().GetZ() + primary_steps.first*primary->GetDir().GetZ())));
		rows->Set<float>("cos_theta", float(cos(primary->GetDir().GetZenith())));
	}
	
	uint32_t m = 0;
	float *energies = rows->GetPointer<float>("energy");
	float *radii = rows->GetPointer<float>("radius");
	
	log_trace("%zu total tracks", Track::Harvest(mctree, *mmctracks).size());
	
	BOOST_FOREACH(const Track &track, Track::Harvest(mctree, *mmctracks)) {
		// MuonGun bundles are attached directly to the primary
		if (mctree.depth(track) > 1)
			continue;
		std::pair<double, double> steps =
		    surface_->GetIntersection(track.GetPos(), track.GetDir());
		float energy = float(track.GetEnergy(steps.first));
		log_trace("energy after %f m: %.1e", steps.first, energy);
		if (energy > 0) {
			if (m < maxMultiplicity_) {
				energies[m] = energy;
				radii[m] = float(GetRadius(*primary, track.GetPos(steps.first)));
			}
			m++;
		}
	}
	
	rows->Set("multiplicity", m);
	
	return 1;
}

/**
 * @brief Interface between WeightCalculator and IceTray
 *
 * WeightCalculatorModule handles the details of extracting energies and
 * radial offsets of muons from an I3MCTree and MMCTrackList.
 */
class WeightCalculatorModule : public I3Module, protected WeightCalculator {
public:
	WeightCalculatorModule(const I3Context &ctx) : I3Module(ctx)
	{
		AddOutBox("OutBox");
		AddParameter("Model", "Muon flux model for which to calculate a weight", boost::shared_ptr<BundleModel>());
		AddParameter("Generator", "Generation spectrum for the bundles to be weighted", generator_);
	}
	
	void Configure()
	{
		boost::shared_ptr<BundleModel> model;
		GetParameter("Model", model);
		GetParameter("Generator", generator_);
		
		if (!model)
			log_fatal("No flux model configured!");
		flux_ = model->flux;
		radius_ = model->radius;
		energy_ = model->energy;
		
		if (!generator_)
			log_fatal("No generator configured!");
		
		surface_ = generator_->GetInjectionSurface();
		if (!surface_)
			log_fatal("No surface configured!");
	}
	
	void DAQ(I3FramePtr frame)
	{
		// First, harvest the muons in the bundle at their points of injection, storing
		// everything that's necessary to estimate the energy lost up to an arbitrary point
		I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
		I3MMCTrackListConstPtr mmctracks = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");
		if (!mctree)
			log_fatal("I3MCTree missing!");
		// if (!mmctracks)
		// 	log_fatal("I3MMCTrackList missing!");
		
		const I3MCTree::const_iterator primary = mctree->begin();
		std::pair<double, double> steps =
		    surface_->GetIntersection(primary->GetPos(), primary->GetDir());
		BundleConfiguration bundlespec;
		
		if (mmctracks) {
			std::list<Track> tracks = Track::Harvest(*mctree, *mmctracks);
			BOOST_FOREACH(const Track &track, tracks) {
				// Omit secondary muons
				boost::optional<I3Particle> parent = mctree->parent(track);
				if (parent && parent->GetType() == I3Particle::NuclInt)
					continue;
				bundlespec.push_back(BundleEntry(
				    GetRadius(*primary, track.GetPos(steps.first)), track.GetEnergy(steps.first)));
			}
		} else {
			// log_warn("No MMCTrackList found in the frame! Assuming that everything starts on the sampling surface...");
			BOOST_FOREACH(const I3Particle &track, std::make_pair(mctree->begin(), mctree->end())) {
				if (track.GetType() == I3Particle::MuMinus || track.GetType() == I3Particle::MuPlus)
					bundlespec.push_back(BundleEntry(
					    GetRadius(*primary, track.GetPos()), track.GetEnergy()));
			}
		}
		
		frame->Put(GetName(), boost::make_shared<I3Double>(GetWeight(*primary, bundlespec)));
		PushFrame(frame);
	}
	
	void Finish();
};

// Out-of-line virtual method definition to force the vtable into this translation unit
void WeightCalculatorModule::Finish() {}

}

I3_MODULE(I3MuonGun::WeightCalculatorModule);
