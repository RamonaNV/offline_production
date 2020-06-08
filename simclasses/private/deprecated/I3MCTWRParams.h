/**
 *
 * Definition of I3MCTWRParams class
 *
 *
 * copyright  (C) 2006
 * the IceCube collaboration
 */

#ifndef I3TWRMCPARAMS_H_INCLUDED
#define I3TWRMCPARAMS_H_INCLUDED

#include "dataclasses/Utility.h"
#include "icetray/OMKey.h"

using namespace std;

static const unsigned i3mctwrparams_version_ = 1;

struct I3MCTWRParams
{

  int stop_delay;
  int DMADD_thresh;
  int TWR_thresh;
  double rel_sens;
  int wf_type;
  double afterpulse_prob;
  double afterpulse_time;
  double noise_rate;
  double amplitude;
  double cable_delay;
  bool optical;
  double peArea;

  I3MCTWRParams() 
  {
    stop_delay=INT_MIN;
    DMADD_thresh=INT_MIN;
    TWR_thresh=INT_MIN;
    rel_sens=NAN;
    wf_type=INT_MIN;
    afterpulse_prob=NAN;
    afterpulse_time=NAN;
    noise_rate=NAN;
    amplitude=NAN;
    cable_delay=NAN;
    optical=false;
    peArea=NAN;
  }

  virtual ~I3MCTWRParams();

  template <class Archive> void serialize(Archive & ar, unsigned version);
};

typedef std::map<OMKey, I3MCTWRParams> I3MCTWRParamsMap;


I3_CLASS_VERSION(I3MCTWRParams, i3mctwrparams_version_);

I3_POINTER_TYPEDEFS(I3MCTWRParams);
I3_POINTER_TYPEDEFS(I3MCTWRParamsMap);

#endif

