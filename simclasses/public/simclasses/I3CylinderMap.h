#ifndef I3CYLINDERMAP_H_INCLUDED 
#define I3CYLINDERMAP_H_INCLUDED


#include <icetray/serialization.h>
#include <dataclasses/Utility.h>
#include <icetray/I3Logging.h>
#include <icetray/I3FrameObject.h>
#include <map>
#include <icetray/OMKey.h>
#include <simclasses/I3ExtraGeometryItemCylinder.h>
#include <dataclasses/I3Map.h>

//typedef std::map<OMKey, I3ExtraGeometryItemCylinder > I3CylinderMap;

typedef I3Map<OMKey, I3ExtraGeometryItemCylinder > I3CylinderMap;

I3_POINTER_TYPEDEFS(I3CylinderMap);

#endif
