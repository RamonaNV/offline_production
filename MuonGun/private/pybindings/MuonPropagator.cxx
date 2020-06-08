/** $Id: MuonPropagator.cxx 100755 2013-03-13 02:07:18Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 100755 $
 * $Date: 2013-03-12 20:07:18 -0600 (Tue, 12 Mar 2013) $
 */

#include "MuonGun/MuonPropagator.h"

void
register_MuonPropagator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<MuonPropagator, boost::shared_ptr<MuonPropagator> >("MuonPropagator",
	    init<const std::string&, double, double, double>((
	    arg("medium")="ice", arg("ecut")=-1, arg("vcut")=-1, arg("rho")=1.0)))
	    .def("propagate", &MuonPropagator::propagate, (arg("particle"), arg("distance"),
	        arg("losses")=boost::shared_ptr<std::vector<I3Particle> >()))
	    .def("stochastic_rate", &MuonPropagator::GetStochasticRate, (arg("energy"), arg("fraction"), arg("type")=I3Particle::MuMinus))
	    .def("total_stochastic_rate", &MuonPropagator::GetTotalStochasticRate, (arg("energy"), arg("type")=I3Particle::MuMinus))
	    .def("set_seed", &MuonPropagator::SetSeed)
	    .staticmethod("set_seed")
	;
	
	class_<Crust, boost::shared_ptr<Crust> >("Crust",
	    init<boost::shared_ptr<MuonPropagator> >())
	    .def("add_layer", &Crust::AddLayer)
	    .def("ingest", &Crust::Ingest)
	;
}
