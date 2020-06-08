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
 * $Id: I3CLSimSimpleGeometryFromI3Geometry.cxx 178169 2019-12-19 21:40:42Z jvansanten $
 *
 * @file I3CLSimSimpleGeometryFromI3Geometry.cxx
 * @version $Revision: 178169 $
 * @date $Date: 2019-12-19 14:40:42 -0700 (Thu, 19 Dec 2019) $
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <icetray/I3Units.h>
#include <icetray/I3Logging.h>
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/geometry/I3OMGeo.h"
#include "dataclasses/geometry/I3ModuleGeo.h"
#include "simclasses/I3ExtraGeometryItemCylinder.h"

#include <stdexcept>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

namespace clsim { namespace defaults {

const std::set<int> I3CLSimSimpleGeometryFromI3Geometry::ignoreStrings;
const std::set<unsigned int> I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDs;
const std::set<std::string> I3CLSimSimpleGeometryFromI3Geometry::ignoreSubdetectors;
const int32_t I3CLSimSimpleGeometryFromI3Geometry::ignoreStringIDsSmallerThan = 1;
const int32_t I3CLSimSimpleGeometryFromI3Geometry::ignoreStringIDsLargerThan = std::numeric_limits<int32_t>::max();
const uint32_t I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDsSmallerThan = 1;
const uint32_t I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDsLargerThan = 60;
const bool I3CLSimSimpleGeometryFromI3Geometry::splitIntoPartsAccordingToPosition=false;
const bool I3CLSimSimpleGeometryFromI3Geometry::useHardcodedDeepCoreSubdetector=false;

}}

I3CLSimSimpleGeometry
I3CLSimSimpleGeometryFromI3Geometry(double OMRadius,
                                    double oversizeFactor,
                                    const I3FramePtr &frame,
                                    const std::set<int> &ignoreStrings,
                                    const std::set<unsigned int> &ignoreDomIDs,
                                    const std::set<std::string> &ignoreSubdetectors,
                                    int32_t ignoreStringIDsSmallerThan,
                                    int32_t ignoreStringIDsLargerThan,
                                    uint32_t ignoreDomIDsSmallerThan,
                                    uint32_t ignoreDomIDsLargerThan,
                                    // keep here for backwards compatibility - it's unused, so don't warn
                                    __attribute__((__unused__)) bool splitIntoPartsAccordingToPosition,
                                    bool useHardcodedDeepCoreSubdetector)
{
    I3CLSimSimpleGeometry simpleGeo(OMRadius*oversizeFactor);

    if (!frame) throw std::runtime_error("Received NULL frame pointer!");
    
    log_debug("Ignoring StringNum<%" PRIi32 ", StringNum>%" PRIi32 ", OMNum<%" PRIu32 ", OMNum>%" PRIu32 ".",
              ignoreStringIDsSmallerThan, ignoreStringIDsLargerThan,
              ignoreDomIDsSmallerThan, ignoreDomIDsLargerThan);

    I3ModuleGeoMapConstPtr moduleGeoMap = frame->Get<I3ModuleGeoMapConstPtr>("I3ModuleGeoMap");
    
    if (!moduleGeoMap) {
        if (frame->Has("I3Geometry")) {
            log_fatal("No I3ModuleGeoMap found in frame. There *is* an I3Geometry object. Please run the \"I3GeometryDecomposer\" module before this.");
        } else {
            log_fatal("No I3ModuleGeoMap found in frame. There does not seem to be a geometry.");
        }
    }

    I3MapModuleKeyStringConstPtr subdetectors = frame->Get<I3MapModuleKeyStringConstPtr>("Subdetectors");
    if (!subdetectors) log_error("No subdetector configuration in frame. Missing a \"Subdetectors\" object. Assuming all modules are on the same detector.");

    typedef I3Map<ModuleKey, I3ExtraGeometryItemCylinder> I3MapModuleKeyCylinder;
    auto cableShadow = frame->Get<boost::shared_ptr<const I3MapModuleKeyCylinder>>("CableShadow");

    double cableRadius = 0;
    BOOST_FOREACH(const I3ModuleGeoMap::value_type &i, *moduleGeoMap)
    {
        const ModuleKey &key = i.first;
        const I3ModuleGeo &geo = i.second;
        
        std::string subdetectorName = "Unknown"; // use this if there is no Subdetectors object
        if (subdetectors) {
            I3MapModuleKeyString::const_iterator subdetector_it =
            subdetectors->find(key);
            if (subdetector_it == subdetectors->end()) {
                log_fatal("ModuleKey(%i/%u) not found in \"Subdetectors\".",
                          key.GetString(), key.GetOM());
            }
            subdetectorName = subdetector_it->second;
        }
        
        int32_t string=key.GetString();
        uint32_t dom=key.GetOM();

        if (useHardcodedDeepCoreSubdetector) {
            // special hack for DeepCore
            if ((subdetectorName=="IceCube") || (subdetectorName=="DeepCore"))
            {
                if ((string>=79) && (string<=86)) // these are the DeepCore strings
                {
                    if (geo.GetPos().GetZ()>-30.*I3Units::m) // z=30m is about halfway between the upper and lower parts of DeepCore
                        subdetectorName="DeepCoreUpper";
                    else
                        subdetectorName="DeepCoreLower";
                }
                else if (string > 86)
                {
                    subdetectorName="PINGU";
                }
            }
        }
        
        if ((string < ignoreStringIDsSmallerThan) ||
            (string > ignoreStringIDsLargerThan) ||
            (dom < ignoreDomIDsSmallerThan) ||
            (dom > ignoreDomIDsLargerThan))
            continue;

        if (ignoreStrings.count(string)!=0) continue;
        if (ignoreDomIDs.count(dom)!=0) continue;
        if (ignoreSubdetectors.count(subdetectorName)!=0) continue;

        
        // sanity check
        if (std::abs(geo.GetRadius()-OMRadius) > 0.001*I3Units::mm)
            log_fatal("This version of clsim does only support DOMs with one single size. Configured size=%fmm, size in geometry=%fmm",
                      OMRadius/I3Units::mm, geo.GetRadius()/I3Units::mm);
        // Every DOM has a cable, even if we don't know where it is. If no 
        // position is given, assume the cable is along +x. 
        double cableAngle = 0;
        if (cableShadow) {
            auto shadow = cableShadow->find(key);
            if (shadow != cableShadow->end()) {
                if (cableRadius == 0) {
                    cableRadius = shadow->second.GetRadius();
                } else if (shadow->second.GetRadius() != cableRadius) {
                    log_fatal_stream("This version of clsim only supports a single cable radius. Configured "<<cableRadius<<" but "<<key<<" has radius "<<shadow->second.GetRadius());
                }
                double gap = shadow->second.GetCenter().Magnitude() - OMRadius - cableRadius;
                if (std::abs(gap) > 1*I3Units::mm) {
                    log_fatal_stream("This version of clsim assumes that the cable runs along the outside of the module, but it is "<<gap/I3Units::mm<<" mm from the surface.");
                }
                cableAngle = shadow->second.GetCenter().GetPhi();
            }
        } 

        simpleGeo.AddModule(string,dom, geo.GetPos().GetX(), geo.GetPos().GetY(), geo.GetPos().GetZ(), subdetectorName, cableAngle);
    }

    if (cableRadius > 0) {
        simpleGeo.SetCableRadius(cableRadius*oversizeFactor);
    }

    return simpleGeo;
}


