/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef I3PARTICLEIDMAP_H_INCLUDED
#define I3PARTICLEIDMAP_H_INCLUDED

#include <dataclasses/physics/I3ParticleID.h>


/**
 * I3ParticleIDMap is used to describe the originating particles of I3MCPulses in
 * an I3MCPulseSeriesMap. For each OMKey a map of I3ParticleID to a list of
 * indices of MCPulses is stored, where the indices refer to the corresponding
 * I3MCPulseSeries within the associated I3MCPulseSeriesMap. By convention, each
 * list of indices is kept in sorted order.
 */
typedef std::map<I3ParticleID, std::vector<uint32_t> > ParticlePulseIndexMap;
typedef I3Map<OMKey,  ParticlePulseIndexMap> I3ParticleIDMap;
I3_POINTER_TYPEDEFS(I3ParticleIDMap);

std::ostream& operator<<(std::ostream&, const ParticlePulseIndexMap&);

#endif
