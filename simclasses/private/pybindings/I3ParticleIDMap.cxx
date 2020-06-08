//
//   Copyright (c) 2013   Alex Olivas
//

#include <simclasses/I3ParticleIDMap.hpp>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_I3ParticleIDMap()
{
  class_<ParticlePulseIndexMap >("ParticlePulseIndexMap")
    .def(dataclass_suite<ParticlePulseIndexMap>())
    ;
  
  class_<I3Map<OMKey, ParticlePulseIndexMap>,
         bases<I3FrameObject>,
         I3ParticleIDMapPtr>("I3ParticleIDMap")
    .def(dataclass_suite<I3Map<OMKey, ParticlePulseIndexMap> >())
    ;
  
  register_pointer_conversions<I3ParticleIDMap>();
}
