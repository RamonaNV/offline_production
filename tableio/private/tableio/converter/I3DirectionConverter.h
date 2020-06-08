/**
 * copyright  (C) 2012
 * The Icecube Collaboration
 *
 * $Id: I3DirectionConverter.h 94948 2012-11-04 16:21:52Z nwhitehorn $
 *
 * @version $Revision: 94948 $
 * @date $LastChangedDate: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
 * @author Martin Wolf <martin.wolf@fysik.su.se> $LastChangedBy: nwhitehorn $
 */

#ifndef TABLEIO_I3DIRECTIONCONVERTER_H_INCLUDED
#define TABLEIO_I3DIRECTIONCONVERTER_H_INCLUDED

#include <string>

#include "tableio/I3Converter.h"
#include "dataclasses/I3Direction.h"

class I3DirectionConverter : public I3ConverterImplementation<I3Direction>
{
    public:
        enum BookRefFrame { car = 0, sph = 1, all = 2 };

        I3DirectionConverter()
          : bookRefFrame_(I3DirectionConverter::car),
            tableColumnNamePrefix_("")
        {}

        //______________________________________________________________________
        /** The bookRefFrame argument specifies the reference frame for which
         *  direction data should be booked.
         *  The tableColumnNamePrefix argument specified the prefix string that
         *  should be prefixed to all table column names.
         */
        I3DirectionConverter(
            BookRefFrame bookRefFrame,
            std::string tableColumnNamePrefix
        )
          : bookRefFrame_(bookRefFrame),
            tableColumnNamePrefix_(tableColumnNamePrefix)
        {}

        //______________________________________________________________________
        I3TableRowDescriptionPtr
        CreateDescription(const I3Direction& direction);

        //______________________________________________________________________
        size_t
        FillRows(const I3Direction& direction, I3TableRowPtr row);

    protected:
        /** The enum value in what reference frame the direction should get
         *  booked.
         */
        BookRefFrame bookRefFrame_;

        /** The prefix string that should be prefixed to all table column names.
         */
        std::string tableColumnNamePrefix_;
};

#endif // TABLEIO_I3DIRECTIONCONVERTER_H_INCLUDED
