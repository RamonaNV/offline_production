/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3WaveformSeriesMapConverter.h 87992 2012-05-06 21:54:21Z jvansanten $
 *
 * @version $Revision: 87992 $
 * @date $LastChangedDate: 2012-05-06 15:54:21 -0600 (Sun, 06 May 2012) $
 * @author Eike Middell <eike.middell@desy.de> $LastChangedBy: jvansanten $
 */

#ifndef TABLEIO_I3WAVEFORMSERIESMAPCONVERTER_H_INCLUDED
#define TABLEIO_I3WAVEFORMSERIESMAPCONVERTER_H_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/physics/I3Waveform.h"

class I3WaveformSeriesMapConverter : public I3ConverterImplementation< I3WaveformSeriesMap > {
public:
    I3WaveformSeriesMapConverter(bool calibrate = false, bool bookGeometry = false);

private:
    I3TableRowDescriptionPtr CreateDescription(const I3WaveformSeriesMap& waveforms);
    size_t FillRows(const I3WaveformSeriesMap& waveforms, I3TableRowPtr rows);
    size_t GetNumberOfRows(const I3WaveformSeriesMap& waveforms);

    bool calibrate_;
    bool bookGeometry_;
};

#endif // TABLEIO_I3WAVEFORMCONVERTER_H_INCLUDED
