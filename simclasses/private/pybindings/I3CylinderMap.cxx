#include <simclasses/I3CylinderMap.h>
#include <icetray/python/dataclass_suite.hpp>

using boost::python::class_;
using boost::python::dataclass_suite;
using boost::python::bases;

void register_I3CylinderMap()
{
  class_<I3CylinderMap,  I3CylinderMapPtr, bases<I3FrameObject> >("I3CylinderMap")
    .def(dataclass_suite<I3CylinderMap >())
    ;
  register_pointer_conversions<I3CylinderMap>();
 
  {
      typedef I3Map<ModuleKey,I3ExtraGeometryItemCylinder> map_t;
      class_<map_t, boost::shared_ptr<map_t>, bases<I3FrameObject> >("I3MapModuleKeyI3ExtraGeometryItemCylinder")
        .def(dataclass_suite<map_t >())
        ;
      register_pointer_conversions<map_t>();
  }

}
