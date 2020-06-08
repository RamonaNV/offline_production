/** $Id: EnergyDistribution.cxx 177771 2019-12-09 09:25:23Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 177771 $
 * $Date: 2019-12-09 02:25:23 -0700 (Mon, 09 Dec 2019) $
 */

#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/RadialDistribution.h>
#include <phys-services/I3RandomService.h>
#include "utils.h"

#include <icetray/python/gil_holder.hpp>

namespace I3MuonGun {

using namespace boost::python;

class PyEnergyDistribution : public EnergyDistribution, public wrapper<EnergyDistribution> {
public:
	virtual double GetLog(double depth, double cos_theta,
	    unsigned multiplicity, double radius, log_value energy) const override
	{
		detail::gil_holder lock;
		return get_override("GetLog")(depth, cos_theta, multiplicity, radius, double(energy));
	}
	virtual std::vector<std::pair<double,double> > Generate(I3RandomService &rng,
	    double depth, double cos_theta,
	    unsigned multiplicity, unsigned nsamples) const override
	{
		detail::gil_holder lock;
		return get_override("Generate")(rng, depth, cos_theta, multiplicity, nsamples);
	}
	virtual bool operator==(const EnergyDistribution&) const override
	{
		return false;
	}
	

};

}

void register_EnergyDistribution()
{
	using namespace boost::python;
	using namespace I3MuonGun;
	
	class_<EnergyDistribution, EnergyDistributionPtr, boost::noncopyable>("EnergyDistribution", no_init)
	    DEF("__call__", &EnergyDistribution::operator(), (arg("depth"), "cos_theta", "multiplicity", "radius", "energy"))
	    .def("generate", (std::vector<std::pair<double,double> > (EnergyDistribution::*)(I3RandomService&,double,double,unsigned,unsigned) const)&EnergyDistribution::Generate, (arg("rng"), arg("depth"), "cos_theta", "multiplicity", "nsamples"))
	    .def("integrate", &EnergyDistribution::Integrate, (arg("depth"), "cos_theta", "multiplicity", "r_min", "r_max", "e_min", "e_max"))
	#define PROPS (Min)(Max)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, EnergyDistribution, PROPS)
	#undef PROPS
		// .def("generate", (std::pair<double,double> (EnergyDistribution::*)(I3RandomService&,const RadialDistribution&,double,double,unsigned) const)&EnergyDistribution::Generate, (arg("rng"), arg("radial_dist"), arg("depth"), "cos_theta", "multiplicity"))
	;
	
	class_<SplineEnergyDistribution, boost::shared_ptr<SplineEnergyDistribution>,
	    bases<EnergyDistribution> >("SplineEnergyDistribution",
	    init<const std::string&, const std::string&>((arg("singles"), "bundles")))
	;
	
	class_<BMSSEnergyDistribution, boost::shared_ptr<BMSSEnergyDistribution>,
	    bases<EnergyDistribution> >("BMSSEnergyDistribution")
	    .def("get_spectrum", &BMSSEnergyDistribution::GetSpectrum, (arg("depth"), arg("cos_theta"), arg("multiplicity"), arg("radius")))
	;
	
	// class_<PyEnergyDistribution, boost::noncopyable>("EnergyDistributionBase")
	//     	    DEF("__call__", &EnergyDistribution::operator(), (arg("depth"), "cos_theta", "multiplicity", "radius", "energy"))
	//     	    .def("generate", &EnergyDistribution::Generate, (arg("rng"), arg("depth"), "cos_theta", "multiplicity", "radius"))
	// 	    ;
	
	class_<OffsetPowerLaw, boost::shared_ptr<OffsetPowerLaw> >("OffsetPowerLaw",
	    init<double,double,double,double>((arg("gamma"), "offset", "min", "max")))
	    DEF("__call__", &OffsetPowerLaw::operator(), (arg("energy")))
	    .def("generate", &OffsetPowerLaw::Generate)
	    DEF("isf", &OffsetPowerLaw::InverseSurvivalFunction, (arg("p")))
	;
}
