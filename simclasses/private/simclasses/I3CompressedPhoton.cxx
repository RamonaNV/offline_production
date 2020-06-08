/**
 * Copyright (c) 2013
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
 * $Id: I3CompressedPhoton.cxx 158928 2017-10-19 23:26:27Z cweaver $
 *
 * @file I3CompressedPhoton.cxx
 * @version $Revision: 158928 $
 * @date $Date: 2017-10-19 17:26:27 -0600 (Thu, 19 Oct 2017) $
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <simclasses/I3CompressedPhoton.h>
#include <boost/foreach.hpp>

I3CompressedPhoton::~I3CompressedPhoton() { }

void I3CompressedPhoton::SetParticleID(const I3Particle& p) { 
    particleID_ = p.GetMinorID(); 
    particleMajorID_ = p.GetMajorID();
}

template <class Archive>
void I3CompressedPhoton::serialize (Archive &ar, unsigned version)
{
    if (version > i3compressedphoton_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CompressedPhoton class.",version,i3compressedphoton_version_);
    
    ar & make_nvp("time", time_);
    ar & make_nvp("weight", weight_);
    ar & make_nvp("wavelength", wavelength_);
    ar & make_nvp("zenith", zenith_);
    ar & make_nvp("azimuth", azimuth_);
    if (version == 0) {
        double x(x_), y(y_), z(z_);
        ar & make_nvp("x", x);
        ar & make_nvp("y", y);
        ar & make_nvp("z", z);
        x_ = x;
        y_ = y;
        z_ = z;
        groupVelocity_ = NAN;
    } else {
        ar & make_nvp("x", x_);
        ar & make_nvp("y", y_);
        ar & make_nvp("z", z_);
        ar & make_nvp("groupVelocity", groupVelocity_);
    }
    ar & make_nvp("particleID", particleID_);
    ar & make_nvp("particleMajorID", particleMajorID_);
}     

I3_SERIALIZABLE(I3CompressedPhoton);
I3_SERIALIZABLE(I3CompressedPhotonSeriesMap);

std::ostream& I3CompressedPhoton::Print(std::ostream& os) const{
  os << "[I3CompressedPhoton: \n"
     << "        Time: " << time_ << '\n'
     << "      Weight: " << weight_ << '\n'
     << "  Wavelength: " << wavelength_ << '\n'
     << "      Zenith: " << zenith_ << '\n'
     << "     Azimuth: " << azimuth_ << '\n'
     << "    Position: (" << x_ << ',' << y_ << ',' << z_ << ")\n"
     << "  Group Vel.: " << groupVelocity_ << '\n'
     << "  ParticleID: " << particleMajorID_ << ',' << particleID_ << ")]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const I3CompressedPhoton& p){
  return(p.Print(os));
}
