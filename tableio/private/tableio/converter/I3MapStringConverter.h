/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3MapStringConverter.h 87992 2012-05-06 21:54:21Z jvansanten $
 *
 * @version $Revision: 87992 $
 * @date $LastChangedDate: 2012-05-06 15:54:21 -0600 (Sun, 06 May 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: jvansanten $
 */

#ifndef TABLEIO_I3MAPSTRINGCONVERTER_H_INCLUDED
#define TABLEIO_I3MAPSTRINGCONVERTER_H_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/I3Map.h"


template<typename ValueType>
class I3MapStringConverter : public I3ConverterImplementation<I3Map<std::string, ValueType> > {

    I3TableRowDescriptionPtr CreateDescription(const I3Map<std::string, ValueType>& map) {
        if (map.size() == 0)
            log_error("Warning:  I3Map<std::string, ...> used for converter configuration is EMPTY, "
                      "nothing will be booked");
        I3TableRowDescriptionPtr desc(new I3TableRowDescription());
        typedef typename I3Map<std::string, ValueType>::const_iterator citer_type;
        for (citer_type citer = map.begin (); citer != map.end (); ++citer) {
            desc->AddField<ValueType>(citer->first, "", "");
        }
        return desc;
    }

    size_t FillRows(const I3Map<std::string, ValueType>& map, I3TableRowPtr rows) {
        typedef typename I3Map<std::string, ValueType>::const_iterator citer_type;
        for (citer_type citer = map.begin (); citer != map.end (); ++citer) {
            rows->Set<ValueType>(citer->first, citer->second);
        }
        return 1;
    }
};


typedef I3MapStringConverter<double>        I3MapStringDoubleConverter;
typedef I3MapStringConverter<int>           I3MapStringIntConverter;
typedef I3MapStringConverter<bool>          I3MapStringBoolConverter;


#endif // TABLEIO_I3MAPSTRINGCONVERTER_H_INCLUDED

