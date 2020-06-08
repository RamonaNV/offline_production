/** $Id: CompactTrack.cxx 128654 2015-02-04 18:34:51Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 128654 $
 * $Date: 2015-02-04 11:34:51 -0700 (Wed, 04 Feb 2015) $
 */

#include <MuonGun/CompactTrack.h>
#include <icetray/serialization.h>
#include <icetray/I3FrameObject.h>

namespace I3MuonGun {

CompactTrack::CompactTrack(const I3Particle &p)
{
	radius_ = hypot(p.GetPos().GetX(), p.GetPos().GetY());
	energy_ = p.GetEnergy();
	time_ = p.GetTime();
	type_ = p.GetType();
}

template <typename Archive>
void
CompactTrack::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Radius", radius_);
	ar & make_nvp("Energy", energy_);
	ar & make_nvp("Time", time_);
	ar & make_nvp("Type", type_);
}

bool
CompactTrack::operator==(const CompactTrack &other) const
{
	return time_==other.time_ && energy_==other.energy_ && type_==other.type_ && radius_==other.radius_;
}


template <typename Archive>
void
TrackBundle::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Map", base_object<std::map<double, std::vector<CompactTrack> > >(*this));
}

// anchor the vtable
TrackBundle::~TrackBundle() {}

}

I3_SERIALIZABLE(I3MuonGun::CompactTrack);
I3_SERIALIZABLE(I3MuonGun::TrackBundle);


