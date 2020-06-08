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
 * $Id: I3CLSimSimpleGeometryTextFile.cxx 178168 2019-12-19 21:40:33Z jvansanten $
 *
 * @file I3CLSimSimpleGeometryTextFile.cxx
 * @version $Revision: 178168 $
 * @date $Date: 2019-12-19 14:40:33 -0700 (Thu, 19 Dec 2019) $
 * @author Claudio Kopper
 */

#include <icetray/I3Logging.h>
#include "clsim/I3CLSimSimpleGeometryTextFile.h"

#include <stdexcept>
#include <fstream>
#include <limits>

#include <boost/lexical_cast.hpp>

namespace clsim { namespace defaults {
const int32_t I3CLSimSimpleGeometryTextFile::ignoreStringIDsSmallerThan = 1;
const int32_t I3CLSimSimpleGeometryTextFile::ignoreStringIDsLargerThan = std::numeric_limits<int32_t>::max();
const uint32_t I3CLSimSimpleGeometryTextFile::ignoreDomIDsSmallerThan = 1;
const uint32_t I3CLSimSimpleGeometryTextFile::ignoreDomIDsLargerThan = 60;
}}

I3CLSimSimpleGeometry
I3CLSimSimpleGeometryTextFile(double OMRadius,
                              const std::string &filename,
                              int32_t ignoreStringIDsSmallerThan,
                              int32_t ignoreStringIDsLargerThan,
                              uint32_t ignoreDomIDsSmallerThan,
                              uint32_t ignoreDomIDsLargerThan
                              )
{
    I3CLSimSimpleGeometry simplegeo(OMRadius);

    std::ifstream inFile;
    inFile.open(filename.c_str(), std::ifstream::in);

    if (inFile.fail()) throw std::runtime_error("Could not open input file");

    int64_t readString; // 1 - 86
    int64_t readDom;    // 1 - 60
    double readx, ready, readz;  // dom x, y, z read from file

    while (inFile >> readString >> readDom >> readx >> ready >> readz)
    {
        int32_t string;
        uint32_t dom;
        double x;
        double y;
        double z;

        try
        {
            string = boost::lexical_cast<int32_t>(readString);
            dom = boost::lexical_cast<uint32_t>(readDom);
            x = boost::lexical_cast<double>(readx);
            y = boost::lexical_cast<double>(ready);
            z = boost::lexical_cast<double>(readz);
        }
        catch(boost::bad_lexical_cast &)
        {
            throw std::runtime_error("Read error (numeric conversion)!");
        }

        if ((string < ignoreStringIDsSmallerThan) ||
            (string > ignoreStringIDsLargerThan) ||
            (dom < ignoreDomIDsSmallerThan) ||
            (dom > ignoreDomIDsLargerThan))
            continue;

        simplegeo.AddModule(string, dom, x, y, z, "default");
    }

    inFile.close();

    return simplegeo;
}
