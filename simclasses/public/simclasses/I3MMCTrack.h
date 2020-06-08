#ifndef I3MMCTRACK_H
#define I3MMCTRACK_H
#include "dataclasses/physics/I3Particle.h"

static const unsigned i3mmctrack_version_ = 1;

class I3MMCTrack 
{
 public:

  I3MMCTrack() :
    Elost(NAN),
    xi(NAN), yi(NAN), zi(NAN), ti(NAN), Ei(NAN),
    xf(NAN), yf(NAN), zf(NAN), tf(NAN), Ef(NAN),
    xc(NAN), yc(NAN), zc(NAN), tc(NAN), Ec(NAN)
    {};

  double GetXi() const { return xi; }
  double GetYi() const { return yi; }
  double GetZi() const { return zi; }
  double GetTi() const { return ti; }
  double GetEi() const { return Ei; }

  double GetXc() const { return xc; }
  double GetYc() const { return yc; }
  double GetZc() const { return zc; }
  double GetTc() const { return tc; }
  double GetEc() const { return Ec; }

  double GetXf() const { return xf; }
  double GetYf() const { return yf; }
  double GetZf() const { return zf; }
  double GetTf() const { return tf; }
  double GetEf() const { return Ef; }

  double GetElost() const { return Elost; }

  void SetEnter(double xi_, double yi_, double zi_, double ti_, double Ei_){
    xi=xi_;
    yi=yi_;
    zi=zi_;
    ti=ti_;
    Ei=Ei_;
  };

  void SetExit(double xf_, double yf_, double zf_, double tf_, double Ef_){
    xf=xf_;
    yf=yf_;
    zf=zf_;
    tf=tf_;
    Ef=Ef_;
  };

  void SetCenter(double xc_, double yc_, double zc_, double tc_, double Ec_){
    xc=xc_;
    yc=yc_;
    zc=zc_;
    tc=tc_;
    Ec=Ec_;
  };

  void SetDepositedEnergy(double Elost_){
    Elost=Elost_;
  };

  const I3Particle& GetI3Particle() const { return particle; };
  I3Particle& GetParticle(){ return particle; };
  void SetParticle(const I3Particle& p){ particle = p; };

  I3Particle particle;
  double Elost;
  // coordinates at entry point [m,m,m,sec,GeV]
  double xi, yi, zi, ti, Ei;

  // coordinates at exit point
  double xf, yf, zf, tf, Ef;

  // coordinates at point of closest approach
  double xc, yc, zc, tc, Ec;

 private:

  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive& ar, unsigned version);

};

std::ostream& operator<<(std::ostream& oss, const I3MMCTrack& t);

I3_POINTER_TYPEDEFS(I3MMCTrack);

typedef I3Vector<I3MMCTrack> I3MMCTrackList;
I3_POINTER_TYPEDEFS(I3MMCTrackList);

I3_CLASS_VERSION(I3MMCTrack,i3mmctrack_version_);
#endif //I3MMCTRACK_H
