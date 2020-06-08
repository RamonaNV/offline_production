/**
 * copyright  (C) 2012
 * the IceCube Collaboration
 * $Id: $
 *
 * @file I3WimpParamsPybindings.h
 * @date $Date: 2012-12-20$
 * @author mzoll <marcel.zoll@fysik.su.se>
 */

#include <string>

#include "icetray/python/stream_to_string.hpp"
#include "dataclasses/ostream_overloads.hpp"
#include "simclasses/I3WimpParams.h"


namespace pyWimpParams {
  // Set the source as SourceType
  void SetSource(I3WimpParams &object, const WimpSim::SourceType source)
    {object.SetSource(source);};

  // Set the source as String
  void SetSourceString(I3WimpParams &object, const std::string source)
    {object.SetSource(source);};

  // Set the source as a uint
  void SetSourceInt(I3WimpParams &object, const uint source)
    {object.SetSource(source);};
  
  // Get the source
  WimpSim::SourceType GetSource(I3WimpParams &object)
    {return object.GetSource();};
    
  // Set the Time from a I3Time
  void SetTime(I3WimpParams &object, const I3Time time)
    {object.SetTime(time);};

  // Set the time as MJD-double
  void SetTimeMJD(I3WimpParams &object, const double mjd)
    {object.SetTime(mjd);};
  
  //Get the Time
  I3Time GetTime(I3WimpParams &object)
    {return object.GetTime();};
}


using namespace boost::python;

void register_I3WimpParams()
{
  class_<I3WimpParams, bases<I3FrameObject>, boost::shared_ptr<I3WimpParams> > wimpparams =
  class_<I3WimpParams, bases<I3FrameObject>, boost::shared_ptr<I3WimpParams> >("I3WimpParams")
    //.def(init<int, double, int, double, double, double, double, double>(args("Source", "Mass", "Channel", "Nu_Weight", "Lep_Weight", "Had_Weight", "Vgen", "Time"),"Constructor for I3WimpParams"))

      
    #define PROPERTIES (Mass)(Channel)(NuWeight)(LepWeight)(HadWeight)(Vgen)
    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, I3WimpParams, PROPERTIES)
    #undef PROPERTIES
    //.def( freeze() )
    
    .add_property("source", &pyWimpParams::GetSource, &pyWimpParams::SetSource)
    .def("set_sourcestring", &pyWimpParams::SetSourceString)
    .def("set_sourceint", &pyWimpParams::SetSourceInt)
    .add_property("time" ,&pyWimpParams::GetTime, &pyWimpParams::SetTime)
    .def("set_timemjd", &pyWimpParams::SetTimeMJD)
    .def("__str__", &stream_to_string<I3WimpParams>)
    ;
  {
    scope wimpparams_scope(wimpparams);

    enum_<WimpSim::SourceType>("SourceType")
      .value("SUN",WimpSim::SUN)
      .value("EARTH",WimpSim::EARTH)
      .value("UNKNOWN",WimpSim::UNKNOWN)
      .export_values()
      ;
      
    enum_<WimpSim::DecayChannel>("DecayChannel")
      .value("unknown",WimpSim::unknown)
      .value("down",WimpSim::down)
      .value("up",WimpSim::up)
      .value("strange",WimpSim::strange)
      .value("charm",WimpSim::charm)
      .value("bottom",WimpSim::bottom)
      .value("top",WimpSim::top)
      .value("gluon",WimpSim::gluon)
      .value("w",WimpSim::w)
      .value("z",WimpSim::z)
      .value("muon",WimpSim::muon)
      .value("tau",WimpSim::tau)
      .value("nue",WimpSim::nue)
      .value("numu",WimpSim::numu)
      .value("nutau",WimpSim::nutau)
      .value("KKDM",WimpSim::KKDM)
      .export_values()
      ;
  }      
  register_pointer_conversions<I3WimpParams>();
}
