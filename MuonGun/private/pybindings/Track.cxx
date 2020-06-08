
#include "MuonGun/Track.h"
#include "simclasses/I3MMCTrack.h"
#include "icetray/python/list_indexing_suite.hpp"
#include "icetray/python/stream_to_string.hpp"

namespace I3MuonGun {

inline bool
operator==(const Track::Checkpoint &p1, const Track::Checkpoint &p2)
{
	return p1.energy == p2.energy && p1.length == p2.length && p1.offset == p2.offset;
}

inline bool
operator==(const Track::LossSum &p1, const Track::LossSum &p2)
{
	return p1.energy == p2.energy && p1.length == p2.length;
}

inline std::ostream&
operator<<(std::ostream &s, const Track::Checkpoint &p1)
{
	s << "Checkpoint(" << p1.length << "," << p1.energy << "," << p1.offset << ")";
	return s;
}

inline std::ostream&
operator<<(std::ostream &s, const Track::LossSum &p1)
{
	s << "LossSum(" << p1.length << "," << p1.energy << ")";
	return s;
}

}

void
register_Track()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	{
		scope track_scope = 
		class_<Track, TrackPtr, bases<I3Particle> >("Track")
		    .def("get_energy", (double (Track::*)(double) const)&Track::GetEnergy)
		    .def("harvest", &Track::Harvest)
		    .staticmethod("harvest")
		    // non-public testing interface
		    .add_property("_checkpoints", &Track::GetCheckpoints)
		    .add_property("_losses", &Track::GetLosses)
		;
	
		class_<Track::Checkpoint>("Checkpoint", no_init)
		    .def_readonly("length", &Track::Checkpoint::length)
		    .def_readonly("energy", &Track::Checkpoint::energy)
		    .def_readonly("offset", &Track::Checkpoint::offset)
		    .def("__repr__", &stream_to_string<Track::Checkpoint>)
		;
	
		class_<std::vector<Track::Checkpoint> >("CheckpointSeries")
		    .def(list_indexing_suite<std::vector<Track::Checkpoint> >())
		;
	
		class_<Track::LossSum>("LossSum", no_init)
		    .def_readonly("length", &Track::LossSum::length)
		    .def_readonly("energy", &Track::LossSum::energy)
		    .def("__repr__", &stream_to_string<Track::LossSum>)
		;
	
		class_<std::vector<Track::LossSum> >("LossSumSeries")
		    .def(list_indexing_suite<std::vector<Track::LossSum> >())
		;
	}
	
	class_<std::list<Track>, boost::shared_ptr<std::list<Track> > >("TrackList")
	    .def(list_indexing_suite<std::list<Track> >())
	;

}
