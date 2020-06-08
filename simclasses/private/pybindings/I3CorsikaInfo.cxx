#include <simclasses/I3CorsikaInfo.h>
//#include "simclasses/converter/I3CorsikaInfoConverter.h"
#include <icetray/python/dataclass_suite.hpp>
#include <tableio/converter/pybindings.h>

namespace bp=boost::python;

void register_I3CorsikaInfo()
{
  boost::python::class_<I3CorsikaInfo, boost::python::bases<I3FrameObject>,
                        boost::shared_ptr<I3CorsikaInfo> >("I3CorsikaInfo")
    .def_readwrite("run_id"         ,&I3CorsikaInfo::run_id)
    .def_readwrite("n_events"       ,&I3CorsikaInfo::n_events)
    .def_readwrite("primary_type"   ,&I3CorsikaInfo::primary_type)
    .def_readwrite("atmosphere"     ,&I3CorsikaInfo::atmosphere)
    .def_readwrite("oversampling"   ,&I3CorsikaInfo::oversampling)
    .def_readwrite("cylinder_height",&I3CorsikaInfo::cylinder_height)
    .def_readwrite("cylinder_radius",&I3CorsikaInfo::cylinder_radius) 
    .def_readwrite("min_zenith"     ,&I3CorsikaInfo::min_zenith)
    .def_readwrite("max_zenith"     ,&I3CorsikaInfo::max_zenith)
    .def_readwrite("min_energy"     ,&I3CorsikaInfo::min_energy)
    .def_readwrite("max_energy"     ,&I3CorsikaInfo::max_energy)
    .def_readwrite("power_law_index",&I3CorsikaInfo::power_law_index)
    .def("__str__", &stream_to_string<I3CorsikaInfo>)
    ;
  
  register_pointer_conversions<I3CorsikaInfo>();

  //I3CONVERTER_NAMESPACE(simclasses);
  //I3CONVERTER_EXPORT_DEFAULT(I3CorsikaInfoConverter, "Dumps I3CorsikaInfo S-Frame object");
}
