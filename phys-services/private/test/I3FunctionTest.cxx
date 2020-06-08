/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: I3FunctionTest.cxx 124369 2014-10-09 16:17:40Z nega $

    @version $Revision: 124369 $
    @date $Date: 2014-10-09 10:17:40 -0600 (Thu, 09 Oct 2014) $
    @author dule

    @todo
*/

#include <I3Test.h>

#include "phys-services/I3Functions.h"
#include <string>
using namespace std;

TEST_GROUP(I3Functions)

TEST(ParseFilename)
{
  string searchpattern(getenv("I3_TESTDATA"));
  searchpattern.append("/amanda/*.f2k");  //4
  searchpattern.append(";");
  searchpattern.append(getenv("I3_TESTDATA"));
  searchpattern.append("/ama-*.*"); // 4

  vector<string> v = I3Functions::ParseFilename(searchpattern);
  for(vector<string>::iterator iter=v.begin(); iter!=v.end(); iter++)
  {
    cout<<*iter<<endl;
  }

  ENSURE(v.size()==8,"Didn't find the expected number of files");
}
