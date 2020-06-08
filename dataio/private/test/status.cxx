/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: status.cxx 169194 2016-06-01 18:27:22Z cweaver $

    @version $Revision: 169194 $
    @date $Date: 2016-06-01 12:27:22 -0600 (Wed, 01 Jun 2016) $

*/

#include <I3Test.h>
#include <fstream>
#include <icetray/serialization.h>

#include <icetray/I3Logging.h>

#include "dataclasses/status/I3DOMStatus.h"
#include "dataclasses/status/I3DetectorStatus.h"
#include "dataclasses/status/I3TriggerStatus.h"


#include "serialization-test.h"

#include <boost/preprocessor.hpp>

using namespace icecube::archive;
using namespace std;

TEST_GROUP(status);

#define TEST_THESE (I3DetectorStatus)

#define NON_I3FO_ITEMS (I3DOMStatus)(I3TriggerStatus)

#define SERIALIZATION_TEST(r,data,t) SERIALIZE(t)

BOOST_PP_SEQ_FOR_EACH(SERIALIZATION_TEST, ~, TEST_THESE);

