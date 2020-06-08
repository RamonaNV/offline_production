/**
   copyright  (C) 2004
   the icecube collaboration
   $Id: retroaction.cxx 169194 2016-06-01 18:27:22Z cweaver $

   @version $Revision: 169194 $
   @date $Date: 2016-06-01 12:27:22 -0600 (Wed, 01 Jun 2016) $

*/

#include <I3Test.h>
#include <fstream>
#include <icetray/serialization.h>

#include <icetray/I3Logging.h>
#include <icetray/Utility.h>

#include "serialization-test.h"
#include "dataio-test.h"

#include <boost/preprocessor.hpp>
#include <boost/foreach.hpp>

using namespace icecube::archive;
using namespace std;

TEST_GROUP(retroaction);

TEST(read)
{
  vector<string> i3files;
  string ports = GetDataDir();

  glob((ports + "/serialization/*/*.i3").c_str(), i3files);
  ENSURE(i3files.size() != 0);
  BOOST_FOREACH(const string& s, i3files)
    {
      log_info("%s", s.c_str());
      I3FramePtr fp = load_i3_file(s);
      ENSURE((bool)fp);
	
      cout << "From " << s << ":\n" << *fp << "\n";
      for (I3Frame::const_iterator iter = fp->begin();
	   iter != fp->end();
	   iter++)
	cout << iter->first << " deserialized\n";
    }
}
