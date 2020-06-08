//
//   Copyright (c) 2018   The IceCube Collaboration 
//   

#include <icetray/python/context_suite.hpp>
#include <sim-services/I3CosmicEventGenerator.h>

namespace bp = boost::python;

void register_I3CosmicEventGenerator()
{
  bp::class_<I3IncrementalEventGeneratorService,
             I3IncrementalEventGeneratorServicePtr,
             boost::noncopyable>
    ("I3IncrementalEventGeneratorService",bp::no_init)
    ;

  bp::class_<NeutrinoSelector,
             NeutrinoSelectorPtr,
             boost::noncopyable>
    ("NeutrinoSelector",bp::no_init)
    ;
  
  bp::class_<I3CosmicEventGenerator, boost::shared_ptr<I3CosmicEventGenerator>, boost::noncopyable>
    ("I3CosmicEventGenerator", "docstring",
     bp::init<I3IncrementalEventGeneratorServicePtr>()
     )
    .def(bp::init<I3IncrementalEventGeneratorServicePtr,
         std::function<bool(const I3Particle &)> >())
    .def(bp::init<I3IncrementalEventGeneratorServicePtr,
         std::function<bool(const I3Particle &)>,
         NeutrinoSelectorPtr>())
    .def(bp::init<I3IncrementalEventGeneratorServicePtr,
         std::function<bool(const I3Particle &)>,
         NeutrinoSelectorPtr,
         I3IncrementalEventGeneratorServicePtr>())    
    .def("Generate",&I3CosmicEventGenerator::Generate)
    ;
}
