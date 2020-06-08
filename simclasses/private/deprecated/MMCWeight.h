/**
    copyright  (C) 2008 
    the icecube collaboration
    $Id: MMCWeight.h 2008-09-04 mdagost
*/

#ifndef MMCWEIGHT_H
#define MMCWEIGHT_H

#include "icetray/I3FrameObject.h"
#include "dataclasses/Utility.h"

/**
 * @brief A class to store weights from MMC when cross-section reweighting is turned on
 */
class MMCWeight : public I3FrameObject
{

 public:
    double weight;
    double distToModIntPoint;

    MMCWeight() :
      weight(NAN),
      distToModIntPoint(NAN)
      {};

    virtual ~MMCWeight() {};

 private:

  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive& ar, unsigned version);

};

I3_POINTER_TYPEDEFS(MMCWeight);

#endif
