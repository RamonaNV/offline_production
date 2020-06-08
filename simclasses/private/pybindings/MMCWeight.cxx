/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id$
 *
 * @version $Revision: 66215 $
 * @date $LastChangedDate: 2010-08-17 18:53:26 +0200 (Tue, 17 Aug 2010) $
 * @author Fabian Kislat <fabian.kislat@desy.de> $LastChangedBy: kislat $
 */

#include <tableio/converter/pybindings.h>
#include "simclasses/MMCWeight.h"
#include "simclasses/converter/MMCWeightConverter.h"

namespace bp = boost::python;

void register_MMCWeight()
{
  bp::class_<MMCWeight, bp::bases<I3FrameObject>, MMCWeightPtr>("MMCWeight")
    .def_readwrite("weight", &MMCWeight::weight)
    .def_readwrite("distToModIntPoint", &MMCWeight::distToModIntPoint)
    ;
  register_pointer_conversions<MMCWeight>();

  I3CONVERTER_NAMESPACE(simclasses);
  I3CONVERTER_EXPORT_DEFAULT(MMCWeightConverter, "Dumps MMCWeight objects");
}
