/**
 * @brief Python bindings for corsika-reader
 *
 * @copyright (C) 2018 The Icecube Collaboration
 * 
 * @file corsika-reader.cxx
 * @author Kevin Meagher
 * @date June 2018
 *
 */

#include "icetray/load_project.h"
#include "corsika-reader/I3CORSIKAReaderUtils.h"
#include "corsika-reader/I3CORSIKAService.h"

namespace bp=boost::python;

BOOST_PYTHON_MODULE(corsika_reader)
{
  bp::import("icecube.icetray");
  bp::import("icecube.sim_services");

  bp::def("CorsikaToPDG",&I3CORSIKAReaderUtils::CorsikaToPDG);
  bp::def("PDGToCorsika",&I3CORSIKAReaderUtils::PDGToCorsika);  
  
#ifdef USE_CORSIKA_CLIENT 
  bp::class_<CorsikaService,
             boost::shared_ptr<CorsikaService>,
             bp::bases<I3IncrementalEventGeneratorService>,
             boost::noncopyable>
    ("CorsikaService",
     bp::init<std::string>()
     )
    .def("StartShower",&CorsikaService::StartShower)
    .def("NextParticle",&CorsikaService::NextParticle)
    .def("EndEvent",&CorsikaService::EndEvent)
    ;
#endif //USE_CORSIKA_CLIENT
}
    
