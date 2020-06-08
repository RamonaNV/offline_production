/** $Id: TrackBinner.cxx 112530 2013-11-01 17:06:01Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 112530 $
 * $Date: 2013-11-01 11:06:01 -0600 (Fri, 01 Nov 2013) $
 */

#include <MuonGun/TrackBinner.h>

using namespace I3MuonGun;
using namespace I3MuonGun::histogram;
namespace bp = boost::python;

static boost::shared_ptr<histogram_base> get_primary(TrackBinner &t) { return t.primary_; }
static boost::shared_ptr<histogram_base> get_multiplicity(TrackBinner &t) { return t.multiplicity_; }
static boost::shared_ptr<histogram_base> get_radius(TrackBinner &t) { return t.radius_; }
static boost::shared_ptr<histogram_base> get_energy(TrackBinner &t) { return t.energy_; }

static boost::shared_ptr<histogram_base>
get_nu(NeutrinoBinner &t, I3Particle::ParticleType pt, unsigned i)
{
	NeutrinoBinner::histmap::iterator it = t.histograms_.find(std::abs(pt));
	if (it == t.histograms_.end() || i >= it->second.size())
		return boost::shared_ptr<histogram_base>();
	return it->second[i];
}

void register_TrackBinner()
{
	bp::class_<TrackBinner>("TrackBinner", bp::init<double,double,unsigned>((
	    bp::arg("mindepth")=1.0, bp::arg("maxdepth")=5.0, bp::arg("steps")=9)))
	    .def("consume", &TrackBinner::Consume)
	    .add_property("primary", &get_primary)
	    .add_property("multiplicity", &get_multiplicity)
	    .add_property("radius", &get_radius)
	    .add_property("energy", &get_energy)
	;
	
	bp::class_<NeutrinoBinner>("NeutrinoBinner")
	    .def("consume", &NeutrinoBinner::Consume)
	    .def("get", &get_nu)
	;
}