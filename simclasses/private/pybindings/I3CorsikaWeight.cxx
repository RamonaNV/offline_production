#include "icetray/python/dataclass_suite.hpp"
#import "simclasses/I3CorsikaWeight.h"

namespace bp=boost::python;

void register_I3CorsikaWeight()
{

  bp::class_<I3CorsikaWeight,
             bp::bases<I3FrameObject>,
             boost::shared_ptr<I3CorsikaWeight> >
    ("I3CorsikaWeight")
    .def_readonly("primary",    &I3CorsikaWeight::primary)
    .def_readonly("bias",       &I3CorsikaWeight::bias)
    .def_readonly("weight",     &I3CorsikaWeight::weight)
    .def_readonly("max_x",      &I3CorsikaWeight::max_x)
    .def(bp::dataclass_suite<I3CorsikaWeight>())
    ;
}

