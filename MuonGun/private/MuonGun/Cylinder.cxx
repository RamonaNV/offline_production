/** $Id: Cylinder.cxx 148884 2016-07-29 06:40:58Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 148884 $
 * $Date: 2016-07-29 00:40:58 -0600 (Fri, 29 Jul 2016) $
 */

#include <MuonGun/Cylinder.h>

namespace I3MuonGun {

std::pair<double, double>
Cylinder::GetZRange() const
{
	return std::make_pair(GetCenter().GetZ() - GetLength()/2., GetCenter().GetZ() + GetLength()/2.);
}

double
Cylinder::GetTopArea() const
{
	return M_PI*GetRadius()*GetRadius();
}

double
Cylinder::GetSideArea() const
{
	return 2*GetRadius()*GetLength();
}

bool
Cylinder::operator==(const SamplingSurface &s) const
{
	const Cylinder *other = dynamic_cast<const Cylinder*>(&s);
	if (!other)
		return false;
	else 
		return (GetRadius() == other->GetRadius() &&
		    GetLength() == other->GetLength() &&
		    GetCenter() == other->GetCenter());
}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 1)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	if (version == 0) {
		// Class hierarchy changed from v0 to v1, so we have to
		// deserialize by hand
		double radius(0), length(0);
		I3Position center;
		ar & make_nvp("SamplingSurface", base_object<I3Surfaces::SamplingSurface>(*this));
		ar & make_nvp("Length", length);
		ar & make_nvp("Radius", radius);
		ar & make_nvp("Center", center);
		SetLength(length);
		SetRadius(radius);
		SetCenter(center);
	} else {
		ar & make_nvp("Base", base_object<Base>(*this));
	}
}

}

// explicitly instantiate the base classes used
template class I3Surfaces::detail::CylinderBase<I3MuonGun::SamplingSurface>;
template class I3MuonGun::detail::UprightSurface<I3Surfaces::detail::CylinderBase<I3MuonGun::SamplingSurface> >;

I3_SERIALIZABLE(I3MuonGun::Cylinder);
