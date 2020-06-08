#ifndef I3MMCTRACKLISTCONVERTER_H_INCLUDED
#define I3MMCTRACKLISTCONVERTER_H_INCLUDED

#include <simclasses/I3MMCTrack.h>
#include <icetray/I3Frame.h>
#include <tableio/converter/I3VectorConverter.h>
#include <tableio/converter/pybindings.h>

struct convert_I3MMCTrack {
    typedef I3MMCTrack booked_type;

    I3ConverterPtr base_;
    I3FramePtr dummyFrame_;

    convert_I3MMCTrack();
    void AddFields(I3TableRowDescriptionPtr desc, const booked_type& track = booked_type());
    void FillSingleRow(const booked_type& track, I3TableRowPtr row);
};

typedef I3VectorConverter< convert_I3MMCTrack > I3MMCTrackListConverter;

#endif  // I3MMCTRACKLISTCONVERTER_H_INCLUDED
