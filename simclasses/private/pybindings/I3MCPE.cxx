//
//   Copyright (c) 2013   Alex Olivas
//

#include <simclasses/I3MCPE.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_I3MCPESeries()
{  
  {
    scope mcpe_scope = 
      class_<I3MCPE, boost::shared_ptr<I3MCPE> >("I3MCPE")
      .def(dataclass_suite<I3MCPE>())
	  .def(init<>())
	  .def(init<uint32_t>())
	  .def(init<uint32_t,double>())
      .def_readwrite("time",&I3MCPE::time)
      .def_readwrite("npe",&I3MCPE::npe)
      .def_readonly("ID",&I3MCPE::ID)
      ;
  }

  class_<I3MCPESeries, I3MCPESeriesPtr>("I3MCPESeries")
    .def(dataclass_suite<I3MCPESeries>())
  ;

  class_<I3MCPESeriesMap, I3MCPESeriesMapPtr, bases<I3FrameObject> >("I3MCPESeriesMap")
    .def(dataclass_suite<I3MCPESeriesMap>())
  ;
  register_pointer_conversions<I3MCPESeriesMap>();
}
