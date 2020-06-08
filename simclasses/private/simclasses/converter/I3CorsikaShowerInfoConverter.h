/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3CorsikaShowerInfoConverter.h 87997 2012-05-07 00:07:52Z jvansanten $
 *
 * @version $Revision: 87997 $
 * @date $LastChangedDate: 2012-05-06 18:07:52 -0600 (Sun, 06 May 2012) $
 * @author Fabian Kislat <fabian.kislat@desy.de> $LastChangedBy: jvansanten $
 */

#ifndef TOPSIMULATOR_I3CORSIKASHOWERINFOCONVERTER_H_INCLUDED
#define TOPSIMULATOR_I3CORSIKASHOWERINFOCONVERTER_H_INCLUDED

#include <tableio/I3Converter.h>
#include "simclasses/I3CorsikaShowerInfo.h"

class I3CorsikaShowerInfoConverter : public I3ConverterImplementation< I3CorsikaShowerInfo > {
public:
    I3CorsikaShowerInfoConverter();

private:
    I3TableRowDescriptionPtr CreateDescription(const I3CorsikaShowerInfo& info);
    size_t FillRows(const I3CorsikaShowerInfo& info, I3TableRowPtr rows);

    size_t nLongSteps_;
};

#endif // TOPSIMULATOR_I3CORSIKASHOWERINFOCONVERTER_H_INCLUDED
