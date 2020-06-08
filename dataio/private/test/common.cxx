/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: common.cxx 169194 2016-06-01 18:27:22Z cweaver $

    @version $Revision: 169194 $
    @date $Date: 2016-06-01 12:27:22 -0600 (Wed, 01 Jun 2016) $

*/

#include <I3Test.h>
#include <fstream>
#include <icetray/serialization.h>

#include <icetray/I3Logging.h>

#include <dataclasses/I3Direction.h>		
#include <icetray/I3Bool.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Time.h>
#include <dataclasses/TriggerKey.h>
#include <icetray/OMKey.h>

#include "serialization-test.h"

#include <boost/preprocessor.hpp>

using namespace icecube::archive;
using namespace std;

TEST_GROUP(common);

#define TEST_THESE \
  (I3Direction) \
  (I3Bool) \
  (I3Double) \
  (I3Position) \
  (I3Time)

#define NON_I3FO_ITEMS \
  (TriggerKey)	       \
  (OMKey)
  
#define SERIALIZATION_TEST(r,data,t) SERIALIZE(t)

BOOST_PP_SEQ_FOR_EACH(SERIALIZATION_TEST, ~, TEST_THESE);

