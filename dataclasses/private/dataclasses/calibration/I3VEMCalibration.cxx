/**
 *
 * @version $Id: I3VEMCalibration.cxx 151966 2016-12-05 19:23:44Z olivas $
 * @file I3VEMCalibration.cxx
 * @date $Date: 2016-12-05 12:23:44 -0700 (Mon, 05 Dec 2016) $
 */


#include <icetray/serialization.h>
#include "dataclasses/calibration/I3VEMCalibration.h"

I3VEMCalibration::~I3VEMCalibration() {}

template <class Archive>
void I3VEMCalibration::serialize(Archive& ar, unsigned version)
{
  if (version>i3vemcalibration_version_)
    log_fatal("Attempting to read version %u from file but running version %u of I3VEMCalibration class.",version,i3vemcalibration_version_);

    ar & make_nvp("pePerVEM",      pePerVEM);
    ar & make_nvp("muPeakWidth",   muPeakWidth);
    ar & make_nvp("hglgCrossOver", hglgCrossOver);
    ar & make_nvp("corrFactor",    corrFactor);
}

std::ostream& operator<<(std::ostream& oss, const I3VEMCalibration& vc)
{
  oss << "[ I3VEMCalibration :: " << std::endl
      << "                pePerVEM  : " << vc.pePerVEM << std::endl
      << "                muPeakWidth  : " << vc.muPeakWidth << std::endl
      << "                hglgCrossOver  : " << vc.hglgCrossOver << std::endl
      << "                corrFactor  : " << vc.corrFactor << std::endl
      << "]" ;
  return oss;
}


I3_SERIALIZABLE(I3VEMCalibration);

