/**
 *  $Id: I3MMCTrack.cxx 150498 2016-10-04 13:05:26Z olivas $
 *  
 *  Copyright (C) 2007, 2008
 *  The IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */
#include <boost/foreach.hpp>

#include <simclasses/I3MMCTrack.h>
#include <icetray/python/dataclass_suite.hpp>

#include <simclasses/converter/I3MMCTrackListConverter.h>
#include <tableio/converter/pybindings.h>
#include <tableio/converter/I3VectorConverter.h>

using namespace boost::python;

std::string pretty_print_list(const I3MMCTrackList& track_list){
  std::stringstream ss;
  ss<<"[";
  BOOST_FOREACH(const I3MMCTrack& t, track_list){
    ss<<t;
  }
  ss<<"]";    
  return ss.str();
}

bool operator==(const I3MMCTrack& lhs, const I3MMCTrack& rhs){
  return 
    lhs.GetXi() == rhs.GetXi() ||
    lhs.GetYi() == rhs.GetYi() ||
    lhs.GetZi() == rhs.GetZi() ||
    lhs.GetTi() == rhs.GetTi() ||
    lhs.GetEi() == rhs.GetEi() ||
    lhs.GetXc() == rhs.GetXc() ||
    lhs.GetYc() == rhs.GetYc() ||
    lhs.GetZc() == rhs.GetZc() ||
    lhs.GetTc() == rhs.GetTc() ||
    lhs.GetEc() == rhs.GetEc() ||
    lhs.GetXf() == rhs.GetXf() ||
    lhs.GetYf() == rhs.GetYf() ||
    lhs.GetZf() == rhs.GetZf() ||
    lhs.GetTf() == rhs.GetTf() ||
    lhs.GetEf() == rhs.GetEf() ;
}

void register_I3MMCTrack()
{
  class_<I3MMCTrack, boost::shared_ptr<I3MMCTrack> >
    ("I3MMCTrack")    
    .def("GetXi", &I3MMCTrack::GetXi )
    .def("GetYi", &I3MMCTrack::GetYi )
    .def("GetZi", &I3MMCTrack::GetZi )
    .def("GetEi", &I3MMCTrack::GetEi )
    .def("GetTi", &I3MMCTrack::GetTi )
    .def("GetXc", &I3MMCTrack::GetXc )
    .def("GetYc", &I3MMCTrack::GetYc )
    .def("GetZc", &I3MMCTrack::GetZc )
    .def("GetEc", &I3MMCTrack::GetEc )
    .def("GetTc", &I3MMCTrack::GetTc )
    .def("GetXf", &I3MMCTrack::GetXf )
    .def("GetYf", &I3MMCTrack::GetYf )
    .def("GetZf", &I3MMCTrack::GetZf )
    .def("GetEf", &I3MMCTrack::GetEf )
    .def("GetTf", &I3MMCTrack::GetTf )
    .def("GetElost", &I3MMCTrack::GetElost )
    .def("GetI3Particle", &I3MMCTrack::GetI3Particle, return_internal_reference<1>() )

    .def("SetParticle", &I3MMCTrack::SetParticle )
    .add_property("particle", make_function(&I3MMCTrack::GetI3Particle, return_internal_reference<1>() ), &I3MMCTrack::SetParticle)
    
    .def_readwrite("xi", &I3MMCTrack::xi)
    .def_readwrite("yi", &I3MMCTrack::yi)
    .def_readwrite("zi", &I3MMCTrack::zi)
    .def_readwrite("ti", &I3MMCTrack::ti)
    .def_readwrite("Ei", &I3MMCTrack::Ei)
    .def_readwrite("xf", &I3MMCTrack::xf)
    .def_readwrite("yf", &I3MMCTrack::yf)
    .def_readwrite("zf", &I3MMCTrack::zf)
    .def_readwrite("tf", &I3MMCTrack::tf)
    .def_readwrite("Ef", &I3MMCTrack::Ef)
    .def_readwrite("xc", &I3MMCTrack::xc)
    .def_readwrite("yc", &I3MMCTrack::yc)
    .def_readwrite("zc", &I3MMCTrack::zc)
    .def_readwrite("tc", &I3MMCTrack::tc)
    .def_readwrite("Ec", &I3MMCTrack::Ec)
    .def_readwrite("Elost", &I3MMCTrack::Elost)
    .def(dataclass_suite<I3MMCTrack>())
    ;

  class_<I3MMCTrackList, bases<I3FrameObject> >("I3MMCTrackList")
    .def(dataclass_suite<I3MMCTrackList >())
    .def("__str__", &pretty_print_list)
    ;

  register_pointer_conversions<I3MMCTrackList>();

  I3CONVERTER_NAMESPACE(simclasses);
  I3CONVERTER_EXPORT_DEFAULT(I3MMCTrackListConverter, "Converts an I3MMCTrackList (I3Vector of I3MMCTrackList)");
}

