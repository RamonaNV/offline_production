/** $Id: Surface.cxx 137873 2015-09-23 13:50:18Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137873 $
 * $Date: 2015-09-23 07:50:18 -0600 (Wed, 23 Sep 2015) $
 */

#include <MuonGun/Cylinder.h>
#include <MuonGun/ExtrudedPolygon.h>

#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <phys-services/I3RandomService.h>
#include <icetray/python/dataclass_suite.hpp>
#include <MuonGun/Flux.h>
#include <icetray/python/gil_holder.hpp>
#include "utils.h"

static double IntegrateFlux(const I3MuonGun::SamplingSurface &s, I3MuonGun::FluxPtr flux, unsigned m, double cosMin, double cosMax)
{
	return s.IntegrateFlux(boost::bind(boost::cref(*flux), _1, _2, m), cosMin, cosMax);
}

void register_Surface()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<SamplingSurface, SamplingSurfacePtr, bases<I3Surfaces::SamplingSurface>, boost::noncopyable>("SamplingSurface", no_init)
	    .def("integrate_flux", &IntegrateFlux, (arg("self"), arg("flux"), arg("m")=1u, arg("cosMin")=0, arg("cosMax")=1))
	    DEF("acceptance", &SamplingSurface::GetAcceptance, (arg("cosMin")=0., arg("cosMax")=1.))
	;
	
	implicitly_convertible<SamplingSurfacePtr, SamplingSurfaceConstPtr>();
	
	class_<Cylinder, CylinderPtr, bases<SamplingSurface> >("Cylinder",
	    init<double, double, I3Position>((arg("length"), arg("radius"), arg("center")=I3Position(0,0,0))))
	    .def(copy_suite<Cylinder>())
	    #define PROPS (Length)(Radius)(Center)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, Cylinder, PROPS)
	    #undef PROPS
	    .def(self == self)
	;

	implicitly_convertible<CylinderPtr, CylinderConstPtr>();

	class_<ExtrudedPolygon, ExtrudedPolygonPtr, bases<SamplingSurface> >("ExtrudedPolygon",
	    init<const std::vector<I3Position> &, double>((arg("points"), arg("padding")=0)))
	    .add_property("x", &ExtrudedPolygon::GetX)
	    .add_property("y", &ExtrudedPolygon::GetY)
	    .add_property("z", &ExtrudedPolygon::GetZ)
	;
}
