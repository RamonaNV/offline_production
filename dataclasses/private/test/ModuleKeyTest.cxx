/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: ModuleKeyTest.cxx 94949 2012-11-04 16:40:30Z nwhitehorn $

    @version $Revision: 94949 $
    @date $Date: 2012-11-04 09:40:30 -0700 (Sun, 04 Nov 2012) $
    @author pretz

    @todo
*/

#include <I3Test.h>

#include "dataclasses/ModuleKey.h"
#include <string>
using std::string;
using std::cout;
using std::endl;

TEST_GROUP(ModuleKeyTest);


TEST(comparison_operator)
{
  ENSURE(ModuleKey(1,1) != ModuleKey(1,2),"different keys are different");
  ENSURE(ModuleKey(1,3) != ModuleKey(2,3),"different keys are different");
  ENSURE(ModuleKey(1,0) < ModuleKey(1,1),"operator< works as expected");
  ENSURE(ModuleKey(1,0) < ModuleKey(2,0),"operator< works as expected");
  ENSURE(ModuleKey(1,1) < ModuleKey(2,2),"operator< works as expected");
  
}
