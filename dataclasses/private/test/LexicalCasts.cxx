/**
    copyright  (C) 2012
    the icecube collaboration
    $Id: LexicalCasts.cxx 94948 2012-11-04 16:21:52Z nwhitehorn $

    @version $Revision: 94948 $
    @date $Date: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
*/

#include <I3Test.h>
#include <icetray/I3Logging.h>
#include <icetray/I3Context.h>
#include <dataclasses/ModuleKey.h>
#include <boost/lexical_cast.hpp>
#include <string>
using namespace std;

TEST_GROUP(LexicalCasts);

TEST(ModuleKey_roundtrip)
{
  // ModuleKey goes round-trip
  ModuleKey k(3,14);
  string s = boost::lexical_cast<string>(k);
  cout << s << endl;
  ModuleKey k2 = boost::lexical_cast<ModuleKey>(s);
  ENSURE_EQUAL(k, k2);

  k = ModuleKey(-13,0);
  s = boost::lexical_cast<string>(k);
  k2 = boost::lexical_cast<ModuleKey>(s);
  ENSURE_EQUAL(k, k2);
}
