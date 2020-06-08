/** $Id: Flux.cxx 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#include <MuonGun/Flux.h>
#include <icetray/python/gil_holder.hpp>
#include "utils.h"

namespace I3MuonGun {

using namespace boost::python;

class PyFlux : public Flux, public wrapper<Flux> {
public:
	virtual double GetLog(double h, double ct, unsigned m) const override
	{
		detail::gil_holder lock;
		return get_override("GetLog")(h, ct, m);
	}
};

}

void register_Flux()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<Flux, boost::shared_ptr<Flux>, boost::noncopyable>("Flux", no_init)
	    DEF("__call__", &Flux::operator(), (args("depth"), "cos_theta", "multiplicity"))
	    #define PROPS (MinMultiplicity)(MaxMultiplicity)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, Flux, PROPS)
	    #undef PROPS
	;
	
	class_<SplineFlux, bases<Flux> >("SplineFlux", init<const std::string&, const std::string&>())
	;
	
	class_<BMSSFlux, bases<Flux> >("BMSSFlux")
	;
	
	class_<PyFlux, boost::noncopyable>("FluxBase")
    	    DEF("__call__", &Flux::operator(), (args("depth"), "cos_theta", "multiplicity"))
    	    #define PROPS (MinMultiplicity)(MaxMultiplicity)
    	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, Flux, PROPS)
    	    #undef PROPS
	;
}
