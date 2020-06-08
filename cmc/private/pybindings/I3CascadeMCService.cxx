//
//   Copyright (c) 2004, 2005, 2006, 2007   Troy D. Straszheim  
//   
//   $Id$
//
//   This file is part of IceTray.
//
//   IceTray is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <cmc/I3CascadeMCService.h>

using namespace boost::python;

void register_I3CascadeMCService()
{
  class_<I3CascadeMCService, boost::shared_ptr<I3CascadeMCService>, 
         bases<I3PropagatorService>, boost::noncopyable>
    ("I3CascadeMCService",init< I3RandomServicePtr > () )
    .def("SetEnergyThresholdMuons", &I3CascadeMCService::SetEnergyThresholdMuons)
    .def("SetMaxMuons", &I3CascadeMCService::SetMaxMuons)
    .def("SetThresholdSplit", &I3CascadeMCService::SetThresholdSplit)
    .def("SetEnergyThresholdSimulation", &I3CascadeMCService::SetEnergyThresholdSimulation)
    .def("SetStepWidth", &I3CascadeMCService::SetStepWidth)
    ;
}
