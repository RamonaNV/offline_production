/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3ShadowedPhotonRemoverModule.cxx 180411 2020-06-03 16:03:01Z olivas $
 *
 * @file I3ShadowedPhotonRemoverModule.cxx
 * @version $Revision: 180411 $
 * @date $Date: 2020-06-03 10:03:01 -0600 (Wed, 03 Jun 2020) $
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>
#include <vector>
#include <iostream>

#include "clsim/shadow/I3ShadowedPhotonRemoverModule.h"
#include "simclasses/I3ExtraGeometryItemCylinder.h"
#include "simclasses/I3CompressedPhoton.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/I3Constants.h"
#include "simclasses/I3CylinderMap.h"

I3_MODULE(I3ShadowedPhotonRemoverModule);

/**
 * Given a cylinder and a photon determine whether the photon intersects the cylinder
 * backpropagating a given distance.
 */
bool I3ShadowPhotonRemover::is_shadowed(const I3Position& photon_position,
					const I3Direction& photon_direction,
					const I3ExtraGeometryItemCylinder& cylinder,
					const double distance){

  // The zenith and azimuth denotes the direction
  // the particle came from, so using these for phi
  // and theta amounts to back-propagating the photon,
  // which is what we want to do.
  double phi{photon_direction.GetAzimuth()};
  double theta{photon_direction.GetZenith()};
  double dx{distance * sin(theta) * cos(phi)};
  double dy{distance * sin(theta) * sin(phi)};
  double dz{distance * cos(theta)};
  const I3Position& line_end{dx, dy, dz};
  return cylinder.DoesLineIntersect(photon_position, line_end);  
}

I3ShadowedPhotonRemoverModule::I3ShadowedPhotonRemoverModule(const I3Context& context) 
: I3ConditionalModule(context)
{
    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputPhotonSeriesMapName_=inputPhotonSeriesMapName_+"Shadowed";
    AddParameter("OutputPhotonSeriesMapName",
                 "Name of the output I3PhotonSeriesMap frame object.",
                 outputPhotonSeriesMapName_);

    AddParameter("CableMapName",
		 "Map containing all the cables found in the geometry",
		 cylinder_map_name_);

    distance_ = 10.0;
    AddParameter("Distance" ,
		 "Distance from where photon hits DOM to extended distance to last scatter" , 
		 distance_) ;
}

I3ShadowedPhotonRemoverModule::~I3ShadowedPhotonRemoverModule()
{
  log_trace("%s", __PRETTY_FUNCTION__);
}


void I3ShadowedPhotonRemoverModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputPhotonSeriesMapName", outputPhotonSeriesMapName_);
    GetParameter("CableMapName", cylinder_map_name_);
    GetParameter("Distance", distance_);
}


void I3ShadowedPhotonRemoverModule::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3CompressedPhotonSeriesMapConstPtr input = frame->Get<I3CompressedPhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if(!input)
      log_fatal("Frame does not contain an I3CompressedPhotonSeriesMap named \"%s\".",
		inputPhotonSeriesMapName_.c_str());

    const I3CylinderMap& cylinder_map = frame->Get<I3CylinderMap>(cylinder_map_name_);
    const I3Geometry& geo = frame->Get<I3Geometry>();
    
    I3CompressedPhotonSeriesMapPtr output(new I3CompressedPhotonSeriesMap());    
    for(const auto& element: *input){
        const ModuleKey& mkey = element.first;
	const OMKey& omkey{mkey.GetString(), mkey.GetOM()};
        const I3CompressedPhotonSeries& photons = element.second;
	const auto& map_pair = cylinder_map.find(omkey);
	
	if(map_pair == end(cylinder_map))
	  continue;
	
	auto cylinder = map_pair->second;
        I3CompressedPhotonSeries unshadowed_photons;
	auto dom_pos = geo.omgeo.at(omkey).position;
        for(auto photon: photons){
	  // Want the photon position in the detector frame
	  // The position of an I3CompressedPhoton is in the DOM frame.
	  auto photon_pos = photon.GetPos();
	  const auto& lab_frame_position{dom_pos + photon_pos};
	  if(!I3ShadowPhotonRemover::is_shadowed(lab_frame_position,
						 photon.GetDir(),
						 cylinder,
						 distance_)){
	    unshadowed_photons.push_back(photon);	    
	  }
        }

	// if there are photons that survive then add them to the output map
	if(unshadowed_photons.size()){	  
	  (*output)[mkey] = unshadowed_photons;
	}
    }

    frame->Put(outputPhotonSeriesMapName_, output);    
    PushFrame(frame);
}
  
