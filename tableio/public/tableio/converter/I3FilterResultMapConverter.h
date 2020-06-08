/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3FilterResultMapConverter.h 136657 2015-08-20 18:50:01Z kkrings $
 *
 * @version $Revision: 136657 $
 * @date $LastChangedDate: 2015-08-20 12:50:01 -0600 (Thu, 20 Aug 2015) $
 * @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: kkrings $
 */

#ifndef TABLEIO_I3FILTERRESULTMAPCONVERTER_H_INCLUDED
#define TABLEIO_I3FILTERRESULTMAPCONVERTER_H_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/physics/I3FilterResult.h"

class I3FilterResultMapConverter : public I3ConverterImplementation<I3FilterResultMap> {
    private:
        I3TableRowDescriptionPtr CreateDescription(const I3FilterResultMap& frmap);
        size_t FillRows(const I3FilterResultMap& frmap, I3TableRowPtr rows);
};

#endif // TABLEIO_I3FILTERRESULTMAPCONVERTER_H_INCLUDED
