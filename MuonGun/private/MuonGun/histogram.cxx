/** $Id: histogram.cxx 128654 2015-02-04 18:34:51Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 128654 $
 * $Date: 2015-02-04 11:34:51 -0700 (Wed, 04 Feb 2015) $
 */

#include <MuonGun/histogram.h>

namespace I3MuonGun {

histogram::binning::scheme::~scheme() {}
histogram::binning::general::~general() {}

histogram::histogram_base::~histogram_base() {}

}
