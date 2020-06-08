/**
 * copyright  (C) 2004
 * the icecube collaboration
 * $Id: I3GCDFileService.h 94949 2012-11-04 16:40:30Z nwhitehorn $
 *
 * @file I3GCDFileService.h
 * @version $Revision: 94949 $
 * @date $Date: 2012-11-04 09:40:30 -0700 (Sun, 04 Nov 2012) $
 * @author pretz
 */

#ifndef I3GCDFILESOURCE_H
#define I3GCDFILESOURCE_H

#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <interfaces/I3GeometryService.h>
#include <interfaces/I3CalibrationService.h>
#include <interfaces/I3DetectorStatusService.h>

/**
 * @brief A I3GeometryOrigin which reads the geometry from a GCD File
 */
class I3GCDFileGeometryService : public I3GeometryService
{
 private:
  std::string filename_;
  I3GeometryConstPtr geo_;
 public:
  I3GCDFileGeometryService(const std::string& icefile) :
  filename_(icefile)     
    {}
  virtual ~I3GCDFileGeometryService(){}
  I3GeometryConstPtr GetGeometry(I3Time time);
};

/**
 * @brief A I3CalibrationOrigin which reads the geometry from a GCD File
 */
class I3GCDFileCalibrationService : public I3CalibrationService
{
 private:
  std::string filename_;
  I3CalibrationConstPtr cal_;
 public:
 I3GCDFileCalibrationService(const std::string& icefile) :
  filename_(icefile) 
  {}
  virtual ~I3GCDFileCalibrationService(){}
  I3CalibrationConstPtr GetCalibration(I3Time time);
};

/**
 * @brief A I3DetectorStatusOrigin which reads the geometry from a GCD File
 */
class I3GCDFileDetectorStatusService : public I3DetectorStatusService
{
 private:
  std::string filename_;
  I3DetectorStatusConstPtr stat_;
 public:
  I3GCDFileDetectorStatusService(const std::string& icefile) :
    filename_(icefile) 
    {}
  virtual ~I3GCDFileDetectorStatusService(){}
  I3DetectorStatusConstPtr GetDetectorStatus(I3Time time);
};

#endif
