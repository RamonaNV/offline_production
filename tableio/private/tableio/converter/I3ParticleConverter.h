/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3ParticleConverter.h 87992 2012-05-06 21:54:21Z jvansanten $
 *
 * @version $Revision: 87992 $
 * @date $LastChangedDate: 2012-05-06 15:54:21 -0600 (Sun, 06 May 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: jvansanten $
 */

#ifndef TABLEIO_I3PARTICLECONVERTER_H_INCLUDED
#define TABLEIO_I3PARTICLECONVERTER_H_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/physics/I3Particle.h"

class I3ParticleConverter : public I3ConverterImplementation<I3Particle> {
private:
    I3TableRowDescriptionPtr CreateDescription(const I3Particle& particle);
    size_t FillRows(const I3Particle& particle, I3TableRowPtr rows);
};

#endif // TABLEIO_I3PARTICLECONVERTER_H_INCLUDED
