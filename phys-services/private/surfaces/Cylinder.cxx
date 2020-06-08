/** $Id: Cylinder.cxx 171111 2019-01-31 15:46:46Z kjmeagher $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 171111 $
 * $Date: 2019-01-31 08:46:46 -0700 (Thu, 31 Jan 2019) $
 */

#include <phys-services/surfaces/Cylinder.h>

namespace I3Surfaces {

Cylinder::~Cylinder() {}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

std::ostream& Cylinder::Print(std::ostream& os) const
{
  os << "Cylinder("<<GetLength() <<", "<< GetRadius() << ", " << GetCenter() << ")";
  return os;
}

}

std::ostream& operator<<(std::ostream& oss, const I3Surfaces::Cylinder& p)
{
  return(p.Print(oss));
}

// explicitly instantiate the base classes used
template class I3Surfaces::detail::CylinderBase<I3Surfaces::SamplingSurface>;

I3_SERIALIZABLE(I3Surfaces::Cylinder);
