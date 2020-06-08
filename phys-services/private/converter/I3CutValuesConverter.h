#ifndef I3CUTVALUESCONVERTER_H_INCLUDED
#define I3CUTVALUESCONVERTER_H_INCLUDED

/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3CutValuesConverter.h 136731 2015-08-21 21:52:17Z nega $
 *
 * @version $Revision: 136731 $
 * @date $LastChangedDate: 2015-08-21 15:52:17 -0600 (Fri, 21 Aug 2015) $
 * @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: nega $
 */

#include <tableio/I3Converter.h>
#include <phys-services/I3CutValues.h>

class I3CutValuesConverter : public I3ConverterImplementation<I3CutValues > {

private:
    I3TableRowDescriptionPtr CreateDescription(const I3CutValues& cv);
    size_t FillRows(const I3CutValues& cv, I3TableRowPtr rows);
    
};

#endif  // I3CUTVALUESCONVERTER_H_INCLUDED
