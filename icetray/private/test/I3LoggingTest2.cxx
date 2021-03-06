/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: I3LoggingTest2.cxx 168504 2013-03-29 13:05:45Z tschmidt $

    @version $Revision: 168504 $
    @date $Date: 2013-03-29 07:05:45 -0600 (Fri, 29 Mar 2013) $
    @author troy d. straszheim <troy@resophonic.com>
*/

#include <I3Test.h>
#include <icetray/I3Logging.h>

#include <string>
using std::string;
using std::cout;
using std::endl;

TEST_GROUP(I3LoggingTest2);
TEST(two)
{
  log_trace("here's a trace message");
  log_debug("here's a debug message");
  log_info("here's an info message");
  log_info("here's an notice message");
  log_warn("here's a warn message");
  log_error("here's an error message");
  try {
    log_fatal("here's a fatal message");
  } catch (std::exception& e) {
    // we should be here
  }
}


