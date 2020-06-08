/** $Id: Generator.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#include <MuonGun/Generator.h>
#include <MuonGun/SamplingSurface.h>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
#include <phys-services/I3RandomService.h>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>


namespace I3MuonGun {

GenerationProbability::~GenerationProbability() {}
Generator::~Generator() {}

template <typename Archive>
void
GenerationProbability::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("NEvents", numEvents_);
}

template <typename Archive>
void
Generator::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("GenerationProbability", base_object<GenerationProbability>(*this));
}

double
GenerationProbability::GetLogGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	return std::log(double(numEvents_)) + GetLogGenerationProbability(axis, bundle);
}

double
GenerationProbability::GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	return numEvents_*std::exp(GetLogGenerationProbability(axis, bundle));
}

GenerationProbabilityCollection::GenerationProbabilityCollection(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2)
{
	push_back(p1);
	push_back(p2);
}

void
GenerationProbabilityCollection::push_back(const GenerationProbabilityPtr &other)
{
	BOOST_FOREACH(value_type &p, *this)
		if (p && p->IsCompatible(other)) {
			p = p->Clone();
			p->SetTotalEvents(p->GetTotalEvents() + other->GetTotalEvents());
			return;
		}
	std::vector<GenerationProbabilityPtr>::push_back(other);
}

bool
GenerationProbabilityCollection::IsCompatible(GenerationProbabilityConstPtr) const
{
	log_fatal("I should never be called.");
}

double
GenerationProbabilityCollection::GetLogGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundle) const
{
	// Collect log probabilities from members
	std::vector<double> values(this->size(), -std::numeric_limits<double>::infinity());
	BOOST_FOREACH(const value_type &p, *this)
		if (p)
			values.push_back(p->GetLogGeneratedEvents(axis, bundle));
	
	// Calculate log(sum(exp)) in a numerically stable way
	double bias = *std::max_element(values.begin(), values.end());
	double prob = 0.;
	BOOST_FOREACH(double v, values)
		prob += std::exp(v-bias);
	
	return bias + std::log(prob);
}

SamplingSurfaceConstPtr
GenerationProbabilityCollection::GetInjectionSurface() const
{
	SamplingSurfaceConstPtr surface;
	BOOST_FOREACH(const value_type &p, *this) {
		if (p) {
			if (!surface)
				surface = p->GetInjectionSurface();
			else if (!(*surface == *(p->GetInjectionSurface())))
				log_fatal("All injection surfaces in a collection must be identical!");
		}
	}
	
	return surface;
}

GenerationProbabilityPtr
GenerationProbabilityCollection::Clone() const
{
	return boost::make_shared<GenerationProbabilityCollection>(*this);
}

template <typename Archive>
void
GenerationProbabilityCollection::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("GenerationProbability", base_object<GenerationProbability>(*this));
	ar & make_nvp("Vector", base_object<std::vector<GenerationProbabilityPtr> >(*this));
}

GenerationProbabilityPtr
operator*=(GenerationProbabilityPtr p, double n)
{
	p->SetTotalEvents(p->GetTotalEvents()*n);
	
	return p;
}

GenerationProbabilityPtr
operator*(GenerationProbabilityPtr p, double n)
{
	GenerationProbabilityPtr pn = p->Clone();
	pn *= n;
	
	return pn;
}

GenerationProbabilityPtr
operator*(double n, GenerationProbabilityPtr p)
{
	return p*n;
}

GenerationProbabilityPtr
operator+(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2)
{

	// If one or both is a collection, merge
	boost::shared_ptr<GenerationProbabilityCollection> c1 =
	    boost::dynamic_pointer_cast<GenerationProbabilityCollection>(p1);
	boost::shared_ptr<GenerationProbabilityCollection> c2 =
	    boost::dynamic_pointer_cast<GenerationProbabilityCollection>(p2);
	if (c1 && c2) {
		c1 = boost::make_shared<GenerationProbabilityCollection>(*c1);
		BOOST_FOREACH(const GenerationProbabilityPtr &other, *c2)
			c1->push_back(other);
		return c1;
	} else if (c1) {
		c1 = boost::make_shared<GenerationProbabilityCollection>(*c1);
		c1->push_back(p2);
		return c1;
	} else if (c2) {
		c2 = boost::make_shared<GenerationProbabilityCollection>(*c2);
		c2->push_back(c1);
		return c2;
	} else if (p1->IsCompatible(p2)) {
		// If the two elements are identical, just scale up.
		GenerationProbabilityPtr p = p1->Clone();
		p->SetTotalEvents(p1->GetTotalEvents() + p2->GetTotalEvents());
		return p;
	} else {
		// If all else fails, start a new collection
		return boost::make_shared<GenerationProbabilityCollection>(p1, p2);
	}
}

I3Particle
Generator::CreateParallelTrack(double radius, double azimuth,
    const I3Surfaces::Surface &surface, const I3Particle &axis)
{
	I3Particle track;
	track.SetLocationType(I3Particle::InIce);
	track.SetType(I3Particle::MuMinus);
	track.SetDir(axis.GetDir());
	track.SetSpeed(I3Constants::c);
	track.SetPos(axis.GetPos());
	track.SetTime(axis.GetTime());
	
	if (radius > 0) {
		// Shift the track parallel to the axis
		I3Position offset(radius, 0, 0);
		offset.RotateY(axis.GetDir().GetZenith());
		offset.RotateZ(azimuth);
		offset += axis.GetPos();
		// Find the distance from the offset position to the sampling
		// surface, and shift so that all tracks have their origin
		// on the surface, but remain in a plane
		double shift = 
		    surface.GetIntersection(offset, track.GetDir()).first
		    -surface.GetIntersection(axis.GetPos(), axis.GetDir()).first;
		if (std::isfinite(shift)) {
			track.SetTime(axis.GetTime() + shift/track.GetSpeed());
			track.SetPos(offset + shift*track.GetDir());
		} else {
			track.SetPos(offset);
		}
	}
	
	return track;
}

/**
 * @brief Interface between Generator and IceTray
 */
class GeneratorModule : public I3Module {
public:
	GeneratorModule(const I3Context &ctx) : I3Module(ctx), maxEvents_(0), numEvents_(0)
	{
		AddOutBox("OutBox");
		AddParameter("Generator", "Muon bundle generator", generator_);
		
		mctreeName_ = "I3MCTree";
	}
	
	void Configure()
	{
		GetParameter("Generator", generator_);
		
		rng_ = context_.Get<I3RandomServicePtr>();
		if (!rng_)
			log_fatal("No RandomService configured!");
		maxEvents_ = size_t(std::floor(generator_->GetTotalEvents()));
		
		firstFrame_ = true;
	}
	
	void DAQ(I3FramePtr frame)
	{
		if (firstFrame_) {
			firstFrame_ = false;
			I3FramePtr sframe = boost::make_shared<I3Frame>('S');
			sframe->Put(GetName(), generator_);
			PushFrame(sframe);
		}
		
		I3MCTreePtr mctree = boost::make_shared<I3MCTree>();
		BundleConfiguration bundlespec;
		
		generator_->Generate(*rng_, *mctree, bundlespec);
		
		frame->Put(mctreeName_, mctree);
		
		PushFrame(frame);
		if (++numEvents_ >= maxEvents_)
			RequestSuspension();
	}
	
	void Finish();
private:
	GeneratorPtr generator_;
	I3RandomServicePtr rng_;
	size_t maxEvents_, numEvents_;
	std::string mctreeName_;
	bool firstFrame_;
};

// Out-of-line virtual method definition to force the vtable into this translation unit
void GeneratorModule::Finish() {}

}

I3_SERIALIZABLE(I3MuonGun::GenerationProbability);
I3_SERIALIZABLE(I3MuonGun::GenerationProbabilityCollection);
I3_SERIALIZABLE(I3MuonGun::Generator);
I3_MODULE(I3MuonGun::GeneratorModule);
