/** $Id: Injectors.cxx 159352 2017-11-07 22:43:41Z juancarlos $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 159352 $
 * $Date: 2017-11-07 15:43:41 -0700 (Tue, 07 Nov 2017) $
 */

#include <MuonGun/StaticSurfaceInjector.h>
#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <MuonGun/Floodlight.h>
#include <MuonGun/NaturalRateInjector.h>
#include <icetray/python/function.hpp>

void
register_CanCan()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<StaticSurfaceInjector, bases<Generator> >("StaticSurfaceInjector")
		.def(init<SamplingSurfacePtr, FluxPtr,
		    boost::shared_ptr<OffsetPowerLaw>, RadialDistributionPtr>())
		#define PROPS (Flux)(RadialDistribution)(EnergyDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, StaticSurfaceInjector, PROPS)
		#undef PROPS
		.add_property("total_rate", &StaticSurfaceInjector::GetTotalRate)
	;
	
	class_<NaturalRateInjector, bases<Generator> >("NaturalRateInjector")
		.def(init<SamplingSurfacePtr, FluxPtr,
		    EnergyDistributionPtr>())
		#define PROPS (Flux)(EnergyDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, NaturalRateInjector, PROPS)
		#undef PROPS
		.add_property("total_rate", &NaturalRateInjector::GetTotalRate)
	;

	class_<SurfaceScalingFunction, SurfaceScalingFunctionPtr, boost::noncopyable>("SurfaceScalingFunction", no_init)
		.def("__call__", &SurfaceScalingFunction::GetSurface)
	;
	
	class_<BasicSurfaceScalingFunction, BasicSurfaceScalingFunctionPtr, bases<SurfaceScalingFunction> >("BasicSurfaceScalingFunction")
		.def("SetSideScaling", &BasicSurfaceScalingFunction::SetSideScaling)
		.def("SetCapScaling", &BasicSurfaceScalingFunction::SetCapScaling)
		.def("SetRadiusBounds", &BasicSurfaceScalingFunction::SetRadiusBounds)
		.def("SetZBounds", &BasicSurfaceScalingFunction::SetZBounds)
		.def("SetCenterMin", &BasicSurfaceScalingFunction::SetCenterMin)
		.def("SetCenterMax", &BasicSurfaceScalingFunction::SetCenterMax)
	;
	class_<ConstantSurfaceScalingFunction, ConstantSurfaceScalingFunctionPtr, bases<SurfaceScalingFunction> >("ConstantSurfaceScalingFunction",
	    init<SamplingSurfacePtr>())
	;

	class_<EnergyDependentSurfaceInjector, bases<StaticSurfaceInjector> >("EnergyDependentSurfaceInjector",
	    init<CylinderPtr, FluxPtr, boost::shared_ptr<OffsetPowerLaw>, RadialDistributionPtr, SurfaceScalingFunctionPtr>
	    ((arg("surface")=CylinderPtr(), arg("flux")=FluxPtr(), arg("energy")=boost::shared_ptr<OffsetPowerLaw>(),
	    arg("radius")=RadialDistributionPtr(), arg("scaling")=boost::make_shared<BasicSurfaceScalingFunction>())))
		.def("total_rate", &EnergyDependentSurfaceInjector::GetTotalRate)
		.def("target_surface", &EnergyDependentSurfaceInjector::GetTargetSurface)
		#define PROPS (Scaling)(Flux)(EnergyDistribution)(RadialDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, EnergyDependentSurfaceInjector, PROPS)
		#undef PROPS
	;

	class_<Floodlight, boost::shared_ptr<Floodlight>, bases<Generator> >("Floodlight",
	    init<SamplingSurfacePtr, boost::shared_ptr<OffsetPowerLaw>, double, double>((
	    arg("surface")=CylinderPtr(),
	    arg("energyGenerator")=boost::shared_ptr<OffsetPowerLaw>(),
	    arg("cosMin")=-1.,arg("cosMax")=1.)))
	;
	
}

