/** $Id: ExtrudedPolygon.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#include <MuonGun/ExtrudedPolygon.h>

namespace I3MuonGun {

ExtrudedPolygon::~ExtrudedPolygon() {}

template <typename Archive>
void
ExtrudedPolygon::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

}

template class I3Surfaces::ExtrudedPolygonBase<I3MuonGun::SamplingSurface>;
template class I3MuonGun::detail::UprightSurface<I3Surfaces::ExtrudedPolygonBase<I3MuonGun::SamplingSurface> >;

I3_SERIALIZABLE(I3MuonGun::ExtrudedPolygon);
