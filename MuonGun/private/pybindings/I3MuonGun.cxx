/** $Id: I3MuonGun.cxx 119085 2014-04-21 20:54:36Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 119085 $
 * $Date: 2014-04-21 14:54:36 -0600 (Mon, 21 Apr 2014) $
 */

#include <MuonGun/I3MuonGun.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>

void
register_I3MuonGun()
{
	using namespace I3MuonGun;
	namespace bp = boost::python;
	
	bp::def("depth", &GetDepth, "Convert a z coordinate to a depth.");
}
