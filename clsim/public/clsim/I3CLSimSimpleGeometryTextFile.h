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
 * $Id: I3CLSimSimpleGeometryTextFile.h 178168 2019-12-19 21:40:33Z jvansanten $
 *
 * @file I3CLSimSimpleGeometryTextFile.h
 * @version $Revision: 178168 $
 * @date $Date: 2019-12-19 14:40:33 -0700 (Thu, 19 Dec 2019) $
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSIMPLEGEOMETRYTEXTFILE_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRYTEXTFILE_H_INCLUDED

#include "clsim/I3CLSimSimpleGeometry.h"

#include <string>

namespace clsim { namespace defaults { 

struct I3CLSimSimpleGeometryTextFile {
    static const int32_t ignoreStringIDsSmallerThan;
    static const int32_t ignoreStringIDsLargerThan;
    static const uint32_t ignoreDomIDsSmallerThan;
    static const uint32_t ignoreDomIDsLargerThan;
};

}}

/**
 * @brief Describes a detector geometry.
 *
 * Reads from a simple text file.
 */
I3CLSimSimpleGeometry I3CLSimSimpleGeometryTextFile(
    double OMRadius, 
    const std::string &filename,
    int32_t ignoreStringIDsSmallerThan=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreStringIDsSmallerThan,
    int32_t ignoreStringIDsLargerThan=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreStringIDsLargerThan,
    uint32_t ignoreDomIDsSmallerThan=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreDomIDsSmallerThan,
    uint32_t ignoreDomIDsLargerThan=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreDomIDsLargerThan
);

#endif //I3CLSIMSIMPLEGEOMETRYTEXTFILE_H_INCLUDED
