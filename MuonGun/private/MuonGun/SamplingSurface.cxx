/** $Id: SamplingSurface.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#include <MuonGun/SamplingSurface.h>
#include <icetray/I3Logging.h>

namespace I3MuonGun {

template <typename Archive>
void
SamplingSurface::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<I3Surfaces::SamplingSurface>(*this));
}

SamplingSurface::~SamplingSurface() {}

}

I3_SERIALIZABLE(I3MuonGun::SamplingSurface);
