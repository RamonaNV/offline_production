/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef I3MCPECONVERTERS_H_INCLUDED
#define I3MCPECONVERTERS_H_INCLUDED

#include <icetray/I3PointerTypedefs.h>
#include <dataclasses/physics/I3MCHit.h>
#include <simclasses/I3MCPE.h>

/**
 * Conversion functions between I3MCHits and I3MCPEs
 * The main thing to keep in mind is that the conversion
 * is asymmetric.
 * 
 * An I3MCPE is only a SPE or RANDOM (in the I3MCHit parlance).
 * So all I3MCPEs can be converted to I3MCHits, but not all
 * I3MCHits can be converted to I3MCPEs.  Namely pre-, late-,
 * and after-pulses are not converted to I3MCPEs.  This is 
 * why HitToPE returns a shared pointer.  If the conversion
 * was not successful a NULL pointer is returned.  Client
 * code needs to check this before dereferencing.
 *
 * For the same reason, in general, the containers won't
 * contain the same number of elements.
 */
namespace PEConversions{
  I3MCHit PEToHit(const I3MCPE&);
  I3MCPEPtr HitToPE(const I3MCHit&);

  void PEToHit(const I3MCPESeries&, I3MCHitSeries&);
  void HitToPE(const I3MCHitSeries&, I3MCPESeries&);

  void PEToHit(const I3MCPESeriesMap&, I3MCHitSeriesMap&);
  void HitToPE(const I3MCHitSeriesMap&, I3MCPESeriesMap&);
  
}

#endif
