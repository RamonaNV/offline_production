/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: calibration.cxx 169194 2016-06-01 18:27:22Z cweaver $

    @version $Revision: 169194 $
    @date $Date: 2016-06-01 12:27:22 -0600 (Wed, 01 Jun 2016) $

*/

#include <I3Test.h>
#include <fstream>
#include <icetray/serialization.h>

#include <icetray/I3Logging.h>

#include <dataclasses/calibration/I3Calibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include <dataclasses/calibration/I3VEMCalibration.h>

#include "serialization-test.h"

#include <boost/preprocessor.hpp>

using namespace icecube::archive;
using namespace std;

TEST_GROUP(calibration);

#define TEST_THESE (I3Calibration)

#define NON_I3FO_ITEMS (I3DOMCalibration)(I3TankCalibration)(I3VEMCalibration)

#define SERIALIZATION_TEST(r,data,t) SERIALIZE(t)

BOOST_PP_SEQ_FOR_EACH(SERIALIZATION_TEST, ~, TEST_THESE);

