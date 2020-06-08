#include <icetray/serialization.h>
#include "simclasses/I3MMCTrack.h"

std::ostream& operator<<(std::ostream& oss, const I3MMCTrack& t){
  oss << "I3MMCTrack = [\n"
      <<" (xi, yi, zi, ti, Ei) = ("
     << t.xi <<" ,"<< t.yi <<" ,"<< t.zi <<" ,"<< t.ti <<" ,"<< t.Ei <<")"<<'\n'
     << " (xc, yc, zc, tc, Ec) = ("
     << t.xc <<" ,"<< t.yc <<" ,"<< t.zc <<" ,"<< t.tc <<" ,"<< t.Ec <<")"<<'\n'
     << " (xf, yf, zf, tf, Ef) = ("
     << t.xf <<" ,"<< t.yf <<" ,"<< t.zf <<" ,"<< t.tf <<" ,"<< t.Ef <<")"<<'\n'
     << " Elost = "<< t.Elost << '\n'
     << " Particle = "<<t.particle<<"]\n";
  return oss;
}


template <class Archive> 
void I3MMCTrack::serialize(Archive& ar, unsigned version)
{

  if (version>i3mmctrack_version_){
    log_fatal("Attempting to read version %u from file but running version %u of I3MMCTrack class.",
	      version,i3particle_version_);
  }

  ar & make_nvp("I3Particle", particle);
  ar & make_nvp("xi", xi);
  ar & make_nvp("yi", yi);
  ar & make_nvp("zi", zi);
  ar & make_nvp("ti", ti);
  ar & make_nvp("Ei", Ei);
  ar & make_nvp("xf", xf);
  ar & make_nvp("yf", yf);
  ar & make_nvp("zf", zf);
  ar & make_nvp("tf", tf);
  ar & make_nvp("Ef", Ef);
  ar & make_nvp("xc", xc);
  ar & make_nvp("yc", yc);
  ar & make_nvp("zc", zc);
  ar & make_nvp("tc", tc);
  ar & make_nvp("Ec", Ec);
  ar & make_nvp("Elost", Elost);

  if(version == 0){
    //these were saved in MMC natural units of seconds
    //this converts to I3Units
    ti *= I3Units::second;
    tf *= I3Units::second;
    tc *= I3Units::second;
  }
}

I3_SERIALIZABLE(I3MMCTrack);
#ifndef OLD_MMCTRACK
I3_SERIALIZABLE(I3MMCTrackList);
#else
typedef I3MMCTrackList MMCTrackList;
I3_SERIALIZABLE(MMCTrackList);
#endif
