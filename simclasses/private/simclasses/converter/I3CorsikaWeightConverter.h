/**
 * copyright  (C) 2020 The Icecube Collaboration
 */

#ifndef I3CORSIKA_WEIGHT_CONVERTER_H_INCLUDED
#define I3CORSIKA_WEIGHT_CONVERTER_H_INCLUDED

#include "tableio/I3Converter.h"
#include "simclasses/I3CorsikaWeight.h"

class I3CorsikaWeightConverter : public I3ConverterImplementation<I3CorsikaWeight> {
private:
    I3TableRowDescriptionPtr CreateDescription(const I3CorsikaWeight& weight);
    size_t FillRows(const I3CorsikaWeight& weight, I3TableRowPtr rows);
};

#endif // I3CORSIKA_WEIGHT_CONVERTER_H_INCLUDED
