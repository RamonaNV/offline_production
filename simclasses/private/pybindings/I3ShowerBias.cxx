/**
 * Copyright (c) 2020 the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * @author Kevin Meagher
 */

#include <icetray/python/dataclass_suite.hpp>
#include <simclasses/I3ShowerBias.h>

namespace bp = boost::python;

void register_I3ShowerBias()
{

  bp::enum_<I3ShowerBias::BiasParticleType>("BiasParticleType")
    .value("Mu",   I3ShowerBias::Mu)
    .value("NuMu", I3ShowerBias::NuMu)
    .value("NuE",  I3ShowerBias::NuE)   
    ;
  
  bp::class_<I3ShowerBias, bp::bases<I3FrameObject>, boost::shared_ptr<I3ShowerBias> >
    ("I3ShowerBias",bp::init<I3ShowerBias::BiasParticleType,double>())
    .def_readwrite("type"  ,&I3ShowerBias::type)
    .def_readwrite("target",&I3ShowerBias::target) 
    .def(bp::dataclass_suite<I3ShowerBias>())
    ;

  bp::class_<I3ShowerBiasMap, bp::bases<I3FrameObject>, I3ShowerBiasMapPtr>("I3ShowerBiasMap")
    .def(bp::dataclass_suite<I3ShowerBiasMap>())
    ;    
  register_pointer_conversions<I3ShowerBias>();
  register_pointer_conversions<I3ShowerBiasMap>();
}
