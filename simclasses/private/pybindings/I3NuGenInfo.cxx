#include <simclasses/I3NuGenInfo.h>
#include <icetray/python/dataclass_suite.hpp>

namespace bp=boost::python;

void register_I3NuGenInfo()
{
  boost::python::class_<I3NuGenInfo, boost::python::bases<I3FrameObject>,
                        boost::shared_ptr<I3NuGenInfo> >("I3NuGenInfo")
    .add_property("run_id"         ,&I3NuGenInfo::run_id)
    .add_property("n_events"       ,&I3NuGenInfo::n_events)
    .add_property("primary_type"   ,&I3NuGenInfo::primary_type)
    .add_property("flavor_fraction",&I3NuGenInfo::flavor_fraction)
    .add_property("cylinder_height",&I3NuGenInfo::cylinder_height)
    .add_property("cylinder_radius",&I3NuGenInfo::cylinder_radius) 
    .add_property("min_zenith"     ,&I3NuGenInfo::min_zenith)
    .add_property("max_zenith"     ,&I3NuGenInfo::max_zenith)
    .add_property("min_azimuth"    ,&I3NuGenInfo::min_azimuth)
    .add_property("max_azimuth"    ,&I3NuGenInfo::max_azimuth)
    .add_property("min_energy"     ,&I3NuGenInfo::min_energy)
    .add_property("max_energy"     ,&I3NuGenInfo::max_energy)
    .add_property("power_law_index",&I3NuGenInfo::power_law_index)
    .def("__str__", &stream_to_string<I3NuGenInfo>)
    ;
  
  register_pointer_conversions<I3NuGenInfo>();
}
