/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3EventHeaderConverter.h 87992 2012-05-06 21:54:21Z jvansanten $
 *
 * @version $Revision: 87992 $
 * @date $LastChangedDate: 2012-05-06 15:54:21 -0600 (Sun, 06 May 2012) $
 * @author Eike Middell <eike.middell@desy.de> $LastChangedBy: jvansanten $
 */

#ifndef TABLEIO_I3EVENTHEADERCONVERTER_HPP_INCLUDED
#define TABLEIO_I3EVENTHEADERCONVERTER_HPP_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/physics/I3EventHeader.h"

class I3EventHeaderConverter : public I3ConverterImplementation<I3EventHeader > {
private:
    I3TableRowDescriptionPtr CreateDescription(const I3EventHeader & params); 
    size_t FillRows(const I3EventHeader& params, I3TableRowPtr rows);
};
    
#endif // TABLEIO_I3EVENTHEADERCONVERTER_HPP_INCLUDED
