/** $Id: Surface.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#include <phys-services/surfaces/Surface.h>

namespace I3Surfaces {

Surface::~Surface() {}

template <typename Archive>
void
Surface::serialize(Archive &ar __attribute__ ((unused)), unsigned version __attribute__ ((unused)))
{}

}

I3_SERIALIZABLE(I3Surfaces::Surface);
