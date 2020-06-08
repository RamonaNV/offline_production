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
 * $Id: I3CLSimSimpleGeometryFromI3Geometry.h 178168 2019-12-19 21:40:33Z jvansanten $
 *
 * @file I3CLSimSimpleGeometryFromI3Geometry.h
 * @version $Revision: 178168 $
 * @date $Date: 2019-12-19 14:40:33 -0700 (Thu, 19 Dec 2019) $
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSIMPLEGEOMETRYFROMI3GEOMETRY_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRYFROMI3GEOMETRY_H_INCLUDED

#include "clsim/I3CLSimSimpleGeometry.h"

#include "icetray/I3Frame.h"

#include <string>
#include <set>

namespace clsim { namespace defaults { 

struct I3CLSimSimpleGeometryFromI3Geometry {
    static const std::set<int> ignoreStrings;
    static const std::set<unsigned int> ignoreDomIDs;
    static const std::set<std::string> ignoreSubdetectors;
    static const int32_t ignoreStringIDsSmallerThan;
    static const int32_t ignoreStringIDsLargerThan;
    static const uint32_t ignoreDomIDsSmallerThan;
    static const uint32_t ignoreDomIDsLargerThan;
    static const bool splitIntoPartsAccordingToPosition;
    static const bool useHardcodedDeepCoreSubdetector;
};

}}

/**
 * @brief Describes a detector geometry.
 *
 * Reads from an I3Geometry frame
 */
I3CLSimSimpleGeometry I3CLSimSimpleGeometryFromI3Geometry(
    double OMRadius, double oversizeFactor,
    const I3FramePtr &frame,
    const std::set<int> &ignoreStrings=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreStrings,
    const std::set<unsigned int> &ignoreDomIDs=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDs,
    const std::set<std::string> &ignoreSubdetectors=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreSubdetectors,
    int32_t ignoreStringIDsSmallerThan=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreStringIDsSmallerThan,
    int32_t ignoreStringIDsLargerThan=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreStringIDsLargerThan,
    uint32_t ignoreDomIDsSmallerThan=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDsSmallerThan,
    uint32_t ignoreDomIDsLargerThan=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDsLargerThan,
    bool splitIntoPartsAccordingToPosition=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::splitIntoPartsAccordingToPosition,
    bool useHardcodedDeepCoreSubdetector=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::useHardcodedDeepCoreSubdetector
);

#endif //I3CLSIMSIMPLEGEOMETRYFROMI3GEOMETRY_H_INCLUDED
