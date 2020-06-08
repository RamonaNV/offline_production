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
 * $Id: I3CLSimSimpleGeometry.h 178169 2019-12-19 21:40:42Z jvansanten $
 *
 * @file I3CLSimSimpleGeometry.h
 * @version $Revision: 178169 $
 * @date $Date: 2019-12-19 14:40:42 -0700 (Thu, 19 Dec 2019) $
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSIMPLEGEOMETRY_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRY_H_INCLUDED

#include <vector>
#include <string>
#include <cstdint>
#include <limits>

#include <icetray/I3PointerTypedefs.h>

/**
 * @brief Describes a detector geometry.
 *
 * The DOM properties have to be explicitly set by the user.
 */

class I3CLSimSimpleGeometry
{
    
public:
    I3CLSimSimpleGeometry(double OMRadius, double CableRadius=0)
        : OMRadius_(OMRadius), CableRadius_(CableRadius)
    {}

    std::size_t size() const {return stringIDs_.size(); }

    /// This is the radius *with* oversizing applied!
    double GetOMRadius() const {return OMRadius_;}
    double GetCableRadius() const {return CableRadius_;}
    void SetOMRadius(double v) {OMRadius_ = v;}
    void SetCableRadius(double v) {CableRadius_ = v;}

    const std::vector<int32_t> &GetStringIDVector() const {return stringIDs_;}
    const std::vector<uint32_t> &GetDomIDVector() const {return domIDs_;}
    const std::vector<double> &GetPosXVector() const {return posX_;}
    const std::vector<double> &GetPosYVector() const {return posY_;}
    const std::vector<double> &GetPosZVector() const {return posZ_;}
    const std::vector<std::string> &GetSubdetectorVector() const {return subdetectors_;}
    const std::vector<double> &GetCableAngleVector() const {return cableAngles_;}

    int32_t GetStringID(std::size_t pos) const {return stringIDs_.at(pos);}
    uint32_t GetDomID(std::size_t pos) const {return domIDs_.at(pos);}
    double GetPosX(std::size_t pos) const {return posX_.at(pos);}
    double GetPosY(std::size_t pos) const {return posY_.at(pos);}
    double GetPosZ(std::size_t pos) const {return posZ_.at(pos);}
    std::string GetSubdetector(std::size_t pos) const {return subdetectors_.at(pos);}
    double GetCableAngle(std::size_t pos) const { return cableAngles_.at(pos);}

    void AddModule(int32_t string, uint32_t om, double x, double y, double z, const std::string &subdetector, double cableAngle=std::numeric_limits<double>::quiet_NaN())
    {
        stringIDs_.push_back(string);
        domIDs_.push_back(om);
        posX_.push_back(x);
        posY_.push_back(y);
        posZ_.push_back(z);
        subdetectors_.push_back(subdetector);
        cableAngles_.push_back(cableAngle);
    }

private:
    double OMRadius_;
    double CableRadius_;

    std::vector<int32_t> stringIDs_;
    std::vector<uint32_t> domIDs_;
    std::vector<double> posX_;
    std::vector<double> posY_;
    std::vector<double> posZ_;
    std::vector<double> cableAngles_;
    std::vector<std::string> subdetectors_;
};

I3_POINTER_TYPEDEFS(I3CLSimSimpleGeometry);

#endif //I3CLSIMSIMPLEGEOMETRY_H_INCLUDED
