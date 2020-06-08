#include "simclasses/converter/I3MMCTrackListConverter.h"

#include <dataclasses/physics/I3Particle.h>

// I3ConverterRegistry::GetDefaultConverter() requires tableio >= r91434

convert_I3MMCTrack::convert_I3MMCTrack() : dummyFrame_(new I3Frame)
{
    base_ = I3ConverterRegistry::GetDefaultConverter<I3Particle>();
    if (!base_)
        log_fatal("No converter registered for I3Particle!");
} 

void convert_I3MMCTrack::AddFields(I3TableRowDescriptionPtr desc, const booked_type& track)
{
    *desc << *(base_->GetDescription(track.GetI3Particle())); // this adds all fields of I3ParticleConverter for the I3Particle
    desc->AddField<double>("Ec", "GeV", "Muon energy at closest point to detector center");
    desc->AddField<double>("Ef", "GeV", "Muon energy at exiting mmc cylinder volume");
    desc->AddField<double>("Ei", "GeV", "Muon energy at entering mmc cylinder volume");
    desc->AddField<double>("Elost", "GeV", "Muon energy loss in mmc cylinder volume");
    desc->AddField<double>("Xc", "m", "Muon x position at closest point to detector center");
    desc->AddField<double>("Xf", "m", "Muon x position at exiting mmc cylinder volume");
    desc->AddField<double>("Xi", "m", "Muon x position at entering mmc cylinder volume");
    desc->AddField<double>("Yc", "m", "Muon y position at closest point to detector center");
    desc->AddField<double>("Yf", "m", "Muon y position at exiting mmc cylinder volume");
    desc->AddField<double>("Yi", "m", "Muon y position at entering mmc cylinder volume");
    desc->AddField<double>("Zc", "m", "Muon z position at closest point to detector center");
    desc->AddField<double>("Zf", "m", "Muon z position at exiting mmc cylinder volume");
    desc->AddField<double>("Zi", "m", "Muon z position at entering mmc cylinder volume");
    desc->AddField<double>("Tc", "ns", "Muon time at closest point to detector center");
    desc->AddField<double>("Tf", "ns", "Muon time at exiting mmc cylinder volume");
    desc->AddField<double>("Ti", "ns", "Muon time at entering mmc cylinder volume");
}

void convert_I3MMCTrack::FillSingleRow(const booked_type& track, I3TableRowPtr row)
{
    size_t nrows = base_->Convert(track.GetI3Particle(), row, dummyFrame_); // this fills all fields of I3ParticleConverter for the I3Particle
    assert(nrows == 1);
    row->Set<double>("Ec", track.GetEc());
    row->Set<double>("Ef", track.GetEf());
    row->Set<double>("Ei", track.GetEi());
    row->Set<double>("Elost", track.GetElost());
    row->Set<double>("Xc", track.GetXc());
    row->Set<double>("Xf", track.GetXf());
    row->Set<double>("Xi", track.GetXi());
    row->Set<double>("Yc", track.GetYc());
    row->Set<double>("Yf", track.GetYf());
    row->Set<double>("Yi", track.GetYi());
    row->Set<double>("Zc", track.GetZc());
    row->Set<double>("Zf", track.GetZf());
    row->Set<double>("Zi", track.GetZi());
    row->Set<double>("Tc", track.GetTc());
    row->Set<double>("Tf", track.GetTf());
    row->Set<double>("Ti", track.GetTi());
} 
