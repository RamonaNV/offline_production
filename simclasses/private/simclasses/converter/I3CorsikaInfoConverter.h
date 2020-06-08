/**
 * copyright  (C) 2020 The Icecube Collaboration
 *
 * $Id: I3CorsikaShowerInfoConverter.h 87997 2012-05-07 00:07:52Z jvansanten $
 *
 * @version $Revision: 87997 $
 * @date $LastChangedDate: 2012-05-06 19:07:52 -0500 (Sun, 06 May 2012) $
 * @author Kevin Meagher $LastChangedBy: jvansanten $
 */

#ifndef I3CORSIKAINFOCONVERTER_H_INCLUDED
#define I3CORSIKAINFOCONVERTER_H_INCLUDED

#include <tableio/I3Converter.h>
#include "simclasses/I3CorsikaInfo.h"

class I3CorsikaInfoConverter : public I3ConverterImplementation< I3CorsikaInfo > {
public:
  //I3CorsikaInfoConverter() : I3ConverterImplementation<I3CorsikaInfo>(){}
private:
  I3TableRowDescriptionPtr CreateDescription(const I3CorsikaInfo& info);
  size_t FillRows(const I3CorsikaInfo& info, I3TableRowPtr rows);
  I3Frame::Stream GetStop() {return I3Frame::Simulation;}
};

#endif // TOPSIMULATOR_I3CORSIKASHOWERINFOCONVERTER_H_INCLUDED
