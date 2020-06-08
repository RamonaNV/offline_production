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
 * $Id: I3ShadowedPhotonRemoverModule.h 180403 2020-06-02 21:36:46Z olivas $
 *
 * @file I3ShadowedPhotonRemoverModule.h
 * @version $Revision: 180403 $
 * @date $Date: 2020-06-02 15:36:46 -0600 (Tue, 02 Jun 2020) $
 * @author Claudio Kopper
 */

#ifndef I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED
#define I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED

#include <string>

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"
#include "dataclasses/geometry/I3Geometry.h"
#include "simclasses/I3CylinderMap.h"

/**
 * @brief This module removes photons that have paths intersecting 
 *   with any shadowing part of the detecor (such as cables).
 *   This code is NOT functional at the moment.
 */
namespace I3ShadowPhotonRemover{
  bool is_shadowed(const I3Position& photon_position,
		   const I3Direction& photon_direction,
		   const I3ExtraGeometryItemCylinder& cylinder,
		   const double distance);
}

class I3ShadowedPhotonRemoverModule : public I3ConditionalModule
{
public:
  I3ShadowedPhotonRemoverModule(const I3Context& ctx);
  ~I3ShadowedPhotonRemoverModule();
  virtual void Configure();
  void DAQ(I3FramePtr frame);
    
private:
  /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
  std::string inputPhotonSeriesMapName_;
  std::string outputPhotonSeriesMapName_;
  std::string cylinder_map_name_;
  double distance_;
  
  I3CylinderMapConstPtr cylinder_map_;
  // default, assignment, and copy constructor declared private
  I3ShadowedPhotonRemoverModule();
  I3ShadowedPhotonRemoverModule(const I3ShadowedPhotonRemoverModule&);
  I3ShadowedPhotonRemoverModule& operator=(const I3ShadowedPhotonRemoverModule&);
  
  SET_LOGGER("I3ShadowedPhotonRemoverModule");
};

#endif //I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED
