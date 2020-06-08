/**
 *  $Id: module.cxx 179712 2020-04-08 15:56:27Z kjmeagher $
 *  
 *  Copyright (C) 2008
 *  Troy D. Straszheim  <troy@icecube.umd.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
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

#include <icetray/load_project.h>
#include <boost/preprocessor.hpp>

#define REGISTER_THESE_THINGS                       \
    (I3Photon)(I3CompressedPhoton)                  \
    (I3Converters)(I3ExtraGeometryItem)		    \

#define I3_REGISTRATION_FN_DECL(r, data, t) void BOOST_PP_CAT(register_,t)();
#define I3_REGISTER(r, data, t) BOOST_PP_CAT(register_,t)();

BOOST_PP_SEQ_FOR_EACH(I3_REGISTRATION_FN_DECL, ~, REGISTER_THESE_THINGS)

void register_I3MMCTrack();
void register_CorsikaLongStep();
void register_I3CorsikaShowerInfo();
void register_I3MCPulse();
void register_I3MCPESeries();
void register_I3WimpParams();
void register_I3ParticleIDMap();
void register_I3CylinderMap();
void register_I3NuGenInfo();
void register_I3CorsikaInfo();
void register_I3CorsikaWeight();
void register_I3ShowerBias();

BOOST_PYTHON_MODULE(simclasses)
{
  load_project("simclasses", false);
  register_I3MMCTrack();
  register_CorsikaLongStep();
  register_I3CorsikaShowerInfo();
  register_I3MCPulse();
  register_I3MCPESeries();
  register_I3WimpParams();
  register_I3ParticleIDMap();
  register_I3CylinderMap();
  register_I3NuGenInfo();  
  register_I3CorsikaInfo();
  register_I3CorsikaWeight();
  register_I3ShowerBias();
  BOOST_PP_SEQ_FOR_EACH(I3_REGISTER, ~, REGISTER_THESE_THINGS);
}

