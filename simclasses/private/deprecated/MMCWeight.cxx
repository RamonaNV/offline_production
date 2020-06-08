#include <icetray/serialization.h>
#include "icetray/I3FrameObject.h"

/**
 * @brief A class to store weights from MMC when cross-section reweighting is turned on
 */
class MMCWeight : public I3FrameObject{

public:
  double weight;
  double distToModIntPoint;
  virtual ~MMCWeight(){};
  std::ostream& Print(std::ostream&) const override;
  
 private:

  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive& ar, unsigned version);

};

template <class Archive>
void MMCWeight::serialize(Archive& ar, unsigned version)
{
  ar & make_nvp("I3FrameObject",     base_object<I3FrameObject>(*this));
  ar & make_nvp("weight",            weight );
  ar & make_nvp("distToModIntPoint", distToModIntPoint );
}

I3_CLASS_VERSION(MMCWeight, 1);
I3_SERIALIZABLE(MMCWeight);

std::ostream& MMCWeight::Print(std::ostream& os) const{
  os << "[MMCWeight:\n"
     << "             Weight: " << weight << '\n'
     << "  DistToModIntPoint: " << distToModIntPoint << "\n]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const MMCWeight& w){
  return(w.Print(os));
}

