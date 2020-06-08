/**
 * copyright  (C) 2004
 * the icecube collaboration
 * $Id: InstallService.cxx 168846 2017-10-18 03:11:03Z olivas $
 *
 * @version $Revision: 168846 $
 * @date $Date: 2017-10-17 21:11:03 -0600 (Tue, 17 Oct 2017) $
 * @author <a href="http://icecube.umd.edu/pretz">John Pretz</a>
 *
 * This is a test which tests that the services are available to modules during
 * the call of the module's destructor
 */

#include <I3Test.h>

#include "icetray/I3Tray.h"
#include "icetray/I3ServiceFactory.h"


TEST_GROUP(InstallService);

// something that installs this service with just one instance for
// every module
struct NameTestService : public I3ServiceFactory
{
  NameTestService(const I3Context& context) 
    : I3ServiceFactory(context) { }
    
  bool InstallService(I3Context& s)
  {
    ENSURE_EQUAL(GetName(), "TESTNAME");
    return true;
  }    
};
  
I3_SERVICE_FACTORY(NameTestService);

TEST(name)
{
  I3Tray tray;
  tray.AddService<NameTestService>("TESTNAME");

  tray.AddModule("BottomlessSource", "bs");
  
  tray.Execute(10);
}

